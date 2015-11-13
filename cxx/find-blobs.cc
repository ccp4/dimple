// g++ -DBUILD_EXE find-blobs.cc -lclipper-ccp4 -lclipper-core -lmmdb2 -o find-blobs

#include "blobs.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <mmdb2/mmdb_manager.h>

using std::string;
using std::vector;
using clipper::Coord_orth;
using clipper::Coord_grid;
using clipper::Coord_frac;
using clipper::RTop_frac;
using clipper::RTop_orth;
using clipper::Xmap;

// popular single-pass algorithm for calculation of variance and mean
class MeanVariance
{
public:
  MeanVariance() : n_(0), mean_(0.), s_(0.) {}

  void add(double x)
  {
    ++n_;
    double delta = x - mean_;
    mean_ += delta / n_;
    s_ += delta * (x - mean_);
  }

  double variance() const { return s_ / (n_ - 1); }
  double stddev() const { return sqrt(variance()); }
  double mean() const { return mean_; }
  int n() const { return n_; }

private:
  int n_;
  double mean_, s_;
};

struct Config
{
  double sigma_level;
  double min_volume;
  double mask_radius; // [A]
  double min_score;
  const char* write_map_to_file;
  const char* write_1_1_map_to_file;
  bool print_mass_center;

  // default values for min_volume and mask_radius
  // are the same as in coot/findligand
  Config() : sigma_level(1.), min_volume(11.), mask_radius(2.), min_score(50),
             write_map_to_file(NULL), write_1_1_map_to_file(NULL),
             print_mass_center(false) {}
};

// unlike Blob, Cluster is for internal use and can contain clipper classes
struct Cluster
{
  // grid points comprising the cluster
  vector<Coord_grid> points;
  // location (.trn() = center, .rot() = ellipsoid shape, not really used now)
  RTop_orth loc;
  MeanVariance mean_variance[3];
  float score; // compatible with coot
};

static
float calculate_cluster_score(const Xmap<float>& xmap, const Cluster& cluster)
{
  float score = 0;
  for (vector<Coord_grid>::const_iterator i =
                  cluster.points.begin(); i != cluster.points.end(); ++i)
    score += xmap.get_data(*i);
  return score;
}

static
bool compare_clusters_by_score(const Cluster &a, const Cluster &b)
{
  return a.score > b.score;
}

static
Coord_orth corth(const Xmap<float>& xmap, const Coord_grid& cg)
{
  return cg.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell());
}

DensityMap::~DensityMap()
{
  delete xmap_;
}

double DensityMap::calculate_stddev() const
{
  return clipper::Map_stats(*xmap_).std_dev();
}

void DensityMap::write_ccp4(const std::string& fname) const
{
  if (xmap_ == NULL)
    return;
  clipper::CCP4MAPfile mapout;
  mapout.open_write(fname);
  mapout.export_xmap(*xmap_);
  mapout.close_write();
}

void DensityMap::mask_by_atoms(mmdb::Manager *mol, double mask_radius)
{
  mmdb::Atom **atoms;
  int n_atoms;
  mol->GetAtomTable(atoms, n_atoms);
  const clipper::Grid_sampling& sampling = xmap_->grid_sampling();

  for (int i = 0; i < n_atoms; i++) {
    const char* name = atoms[i]->residue->name;
    if (strcmp(name, "WAT") == 0 || strcmp(name, "HOH") == 0)
      continue;
    Coord_orth co(atoms[i]->x, atoms[i]->y, atoms[i]->z);
    Coord_frac cf = co.coord_frac(xmap_->cell());
    float ra = mask_radius / xmap_->cell().descr().a();
    float rb = mask_radius / xmap_->cell().descr().b();
    float rc = mask_radius / xmap_->cell().descr().c();
    Coord_frac box0(cf.u() - ra, cf.v() - rb, cf.w() - rc);
    Coord_frac box1(cf.u() + ra, cf.v() + rb, cf.w() + rc);
    clipper::Grid_range gr(box0.coord_grid(sampling),
                           box1.coord_grid(sampling));

    typedef clipper::Xmap_base::Map_reference_coord Mrc;
    for (Mrc iu(*xmap_, gr.min()); iu.coord().u() <= gr.max().u(); iu.next_u())
      for (Mrc iv=iu; iv.coord().v() <= gr.max().v(); iv.next_v())
        for (Mrc iw=iv; iw.coord().w() <= gr.max().w(); iw.next_w()) {
          Coord_orth g = corth(*xmap_, iw.coord());
          if ((g - co).lengthsq() < mask_radius * mask_radius)
            (*xmap_)[iw] = 0.;
        }
  }
}

// all symops x 27 neighbouring cells
static
vector<RTop_frac> alternative_sites(const clipper::Spacegroup& spacegroup)
{
  vector<RTop_frac> rt_ops;
  for (int i = 0; i < spacegroup.num_symops(); ++i) {
    const clipper::Symop& symop = spacegroup.symop(i);
    for (int u = -1; u <= 1; u++)
      for (int v = -1; v <= 1; v++)
        for (int w = -1; w <= 1; w++)
          rt_ops.push_back(RTop_frac(symop.rot(),
                                     symop.trn() + clipper::Vec3<>(u, v, w)));
  }
  return rt_ops;
}

// returns symmetry equivalent location that is closest to the given coords
static
RTop_orth nearest_to_coords(const vector<RTop_orth>& orth_rts,
                            const RTop_orth& loc,
                            const vector<Coord_orth>& coords)
{
  double min_r2 = 1e20;
  vector<RTop_orth>::const_iterator rt;
  for (vector<RTop_orth>::const_iterator i = orth_rts.begin();
                                                  i != orth_rts.end(); ++i) {
    Coord_orth p = Coord_orth(loc.trn()).transform(*i);
    // brute-force search
    for (vector<Coord_orth>::const_iterator j = coords.begin();
                                                     j != coords.end(); ++j) {
      double r2 = (p-*j).lengthsq();
      if (r2 < min_r2) {
        min_r2 = r2;
        rt = i;
      }
    }
  }
  if (min_r2 == 1e20)
    throw std::runtime_error("This should not happen");
  //printf("nearest_to_coords: %s -> %s\n", loc.trn().format().c_str(),
  //                                       (*rt * loc).trn().format().c_str());
  return RTop_orth(*rt * loc);
}

static
void calculate_cluster_props(const Xmap<float>& xmap, Cluster& cluster)
{
  for (vector<Coord_grid>::const_iterator i =
                  cluster.points.begin(); i != cluster.points.end(); ++i) {
    Coord_orth co = corth(xmap, *i);
    cluster.mean_variance[0].add(co.x());
    cluster.mean_variance[1].add(co.y());
    cluster.mean_variance[2].add(co.z());
  }

  clipper::Mat33<double> mat(0,0,0, 0,0,0, 0,0,0);
  Coord_orth ctr(cluster.mean_variance[0].mean(),
                          cluster.mean_variance[1].mean(),
                          cluster.mean_variance[2].mean());
  for (vector<Coord_grid>::const_iterator
          c = cluster.points.begin(); c != cluster.points.end(); ++c) {
    Coord_orth d = corth(xmap, *c) - ctr;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        mat(i,j) += d[i] * d[j];
  }

  // XXX: principal axes of the ellipsoid can be calculated here,
  // but we don't need them yet.

  cluster.loc = RTop_orth(mat, ctr);
}

/*
TODO: Try smarter cluster-finding algorithms,
maybe http://en.wikipedia.org/wiki/OPTICS
*/

// in case of orthorhombic unit cell, returns 18 (=27-1-8) neighbours,
// if angles in the desc arg are far from 90deg the result may be different
static
vector<Coord_grid> get_grid_neighbours(const clipper::Cell_descr& desc)
{
  vector<Coord_grid> nlist;
  nlist.reserve(20);
  clipper::Cell u_cell(clipper::Cell_descr(1., 1., 1.,
                                     desc.alpha(), desc.beta(), desc.gamma()));
  clipper::Grid_sampling u_sampling(1, 1, 1);
  for (int i = -1; i <= 1; ++i)
     for (int j = -1; j <= 1; ++j)
        for (int k = -1; k <= 1; ++k) {
          Coord_grid cg(i,j,k);
          float d2 = cg.coord_frac(u_sampling).lengthsq(u_cell);
          if (d2 > 0.5 && d2 < 2.5)
            nlist.push_back(cg);
        }
  return nlist;
}

vector<Cluster> find_clusters_by_flood_fill(const Xmap<float>& xmap,
                                            float min_volume,
                                            float cut_off)
{
  vector<Cluster> clusters;

  vector<Coord_grid> nabes = get_grid_neighbours(xmap.cell().descr());
  double volume_per_pt = xmap.cell().volume() / xmap.grid_sampling().size();

  // temporary map
  Xmap<float> tmap = xmap;
  for (clipper::Xmap_base::Map_reference_index
                              i = tmap.first(); !i.last(); i.next()) {
    if (tmap[i] > cut_off) {
      Cluster clust;
      clust.points.push_back(i.coord());
      tmap[i] = 0;
      for (size_t j = 0; j < clust.points.size(); ++j) {
        for (size_t k = 0; k < nabes.size(); k++) {
          Coord_grid cg = clust.points[j] + nabes[k];
          if (tmap.get_data(cg) > cut_off) {
            tmap.set_data(cg, 0);
            clust.points.push_back(cg);
          }
        }
      }
      if (clust.points.size() >= 2 &&
          clust.points.size() * volume_per_pt > min_volume)
        clusters.push_back(clust);
    }
  }
  return clusters;
}


mmdb::Manager* read_pdb(const string& pdb_name)
{
  mmdb::InitMatType();
  mmdb::Manager* mm = new mmdb::Manager;
  mm->SetFlag(mmdb::MMDBF_IgnoreBlankLines | mmdb::MMDBF_IgnoreDuplSeqNum |
              mmdb::MMDBF_IgnoreNonCoorPDBErrors | mmdb::MMDBF_IgnoreHash |
              mmdb::MMDBF_IgnoreRemarks);
  int error = mm->ReadCoorFile(pdb_name.c_str());
  if (error)
    std::runtime_error("Error reading PDB file: " + pdb_name);
  mm->PDBCleanup(mmdb::PDBCLEAN_ELEMENT);
  return mm;
}

void DensityMap::read_from_mtz(const string& mtz_name,
                               const string& f_col, const string& phi_col)
{
  clipper::CCP4MTZfile mtzin;
  mtzin.open_read(mtz_name);
  clipper::HKL_info hkl_info;
  // import_hkl_info() copies info from(!) MTZfile
  mtzin.import_hkl_info(hkl_info);
  clipper::MTZcrystal mtz_crystal;
  clipper::HKL_data<clipper::datatypes::F_phi<float> > data(hkl_info,
                                                            mtz_crystal);
  clipper::MTZdataset mtz_dataset;
  // import_hkl_data() only marks data for reading,
  // the actual operation is done in close_read() (!!).
  mtzin.import_hkl_data(data, mtz_dataset, mtz_crystal,
                        "/*/*/["+f_col+" "+phi_col+"]");
  mtzin.close_read();

  clipper::Grid_sampling sampling(hkl_info.spacegroup(), hkl_info.cell(),
                                  hkl_info.resolution(), 1.5);
  delete xmap_;
  xmap_ = new Xmap<float>;
  xmap_->init(hkl_info.spacegroup(), hkl_info.cell(), sampling);
  xmap_->fft_from(data);
}

static
vector<Cluster> find_clusters(const string& pdb_filename,
                              const string& mtz_filename,
                              const Config& config)
{
  string f_col = "FWT";
  string phi_col = "PHWT";

  DensityMap density_map;
  density_map.read_from_mtz(mtz_filename, f_col, phi_col);
  printf("Searching for clusters in density map, using grid: %s\n",
         density_map.xmap()->grid_sampling().format().c_str());
  if (config.write_map_to_file)
    density_map.write_ccp4(config.write_map_to_file);
  float map_stddev = density_map.calculate_stddev();
  if (map_stddev <= 1e-9)
    throw std::runtime_error("Flat density map");

  mmdb::Manager *mm = read_pdb(pdb_filename);

  density_map.mask_by_atoms(mm, config.mask_radius);
  const clipper::Xmap<float>* xmap = density_map.xmap();

  float cut_off = config.sigma_level * map_stddev;
  printf("Density std.dev: %g, cut-off: %g (%g sigma)\n",
         map_stddev, cut_off, config.sigma_level);
  vector<Cluster> all_clusters = find_clusters_by_flood_fill(*xmap,
                                                   config.min_volume, cut_off);

  vector<Cluster> clusters; // for clusters with big enough score
  for (vector<Cluster>::iterator i = all_clusters.begin();
                                 i != all_clusters.end(); ++i) {
    i->score = calculate_cluster_score(*xmap, *i);
    if (i->score >= config.min_score)
      clusters.push_back(*i);
  }

  mmdb::Atom **atoms;
  int n_atoms;
  mm->GetAtomTable(atoms, n_atoms);

  if (config.print_mass_center) {
    mmdb::realtype x, y, z;
    mmdb::GetMassCenter(atoms, n_atoms, x, y, z);
    printf("Protein mass center: %s\n", Coord_orth(x,y,z).format().c_str());
  }

  vector<Coord_orth> pdb_coords;
  pdb_coords.reserve(n_atoms);
  for (int i = 0; i < n_atoms; ++i)
    pdb_coords.push_back(Coord_orth(atoms[i]->x, atoms[i]->y, atoms[i]->z));
  delete mm;

  for (vector<Cluster>::iterator i = clusters.begin(); i != clusters.end(); ++i)
    calculate_cluster_props(*xmap, *i);

  vector<RTop_frac> frac_rts = alternative_sites(xmap->spacegroup());
  vector<RTop_orth> orth_rts(frac_rts.size());
  for (size_t i = 0; i != frac_rts.size(); ++i)
    orth_rts[i] = frac_rts[i].rtop_orth(xmap->cell());
  for (vector<Cluster>::iterator i = clusters.begin(); i != clusters.end(); ++i)
    i->loc = nearest_to_coords(orth_rts, i->loc, pdb_coords);

  std::sort(clusters.begin(), clusters.end(), compare_clusters_by_score);
  return clusters;
}

// for python API
vector<Blob> find_blobs(const string& pdb_filename,
                        const string& mtz_filename,
                        double sigma_level,
                        double min_volume,
                        double mask_radius)
{
  Config config;
  config.sigma_level = sigma_level;
  config.min_volume = min_volume;
  config.mask_radius = mask_radius;
  vector<Cluster> clusters = find_clusters(pdb_filename, mtz_filename, config);
  vector<Blob> blobs(clusters.size());
  for (size_t i = 0; i < clusters.size(); ++i) {
    Blob &b = blobs[i];
    const Cluster &c = clusters[i];
    const clipper::Vec3<>& trn = c.loc.trn();
    b.n_points = c.points.size();
    b.score = c.score;
    b.x = trn[0];
    b.y = trn[1];
    b.z = trn[2];
  }
  return blobs;
}

#if BUILD_EXE

static
const char* usage =
"Usage:  find-blobs [options] file.pdb [file.mtz]\n"
" Options:\n"
"  -s SIGMA        sigma level\n"
"  -w FILE         write 2-1 map (that is used to find blobs) to file\n"
"  -e FILE         write 1-1 map to file\n"
"  -c              print protein mass center\n";

static
bool endswith(std::string const& s, std::string const& p)
{
  return p.size() <= s.size() && std::string(s, s.size() - p.size()) == p;
}

static
int err(const string& e)
{
  fprintf(stderr, "%s\n", e.c_str());
  exit(-1);
}

static
void parse_option(char letter, const char *option, Config* config)
{
  if (letter == 's') {
    char* endptr;
    config->sigma_level = strtod(option, &endptr);
    if (*endptr != '\0')
      err("Error when parsing `" + string(option) + "' as number");
  }
  else if (letter == 'w') {
    config->write_map_to_file = option;
  }
  else if (letter == 'e') {
    config->write_1_1_map_to_file = option;
  }
  else
    err("Unknown option: -" + string(1,letter));
}

int main(int argc, char **argv)
{
  string pdb_filename, mtz_filename;
  Config config;
  bool handle_options = true;
  for (int i = 1; i < argc; ++i) {
    const char *arg = argv[i];
    if (handle_options && arg[0] == '-') {
      // '-' ends options to allow filenames starting with '-'
      if (arg[1] == '\0')
        handle_options = false;
      else if (strcmp(arg+1, "c") == 0) {
        config.print_mass_center = true;
      }
      else if (arg[2] == '\0') {
        if (i+1 == argc)
          err("Unknown option or missing arg for option: " + string(arg));
        parse_option(arg[1], argv[++i], &config);
      }
      else
        parse_option(arg[1], arg+2, &config);
    }
    else if (endswith(arg, ".pdb"))
      pdb_filename = arg;
    else if (endswith(arg, ".mtz"))
      mtz_filename = arg;
    else {
      err("Unknown argument: " + string(arg));
    }
  }
  if (pdb_filename.empty()) {
    fprintf(stderr, "%s", usage);
    return 1;
  }
  if (mtz_filename.empty())
    mtz_filename = pdb_filename.substr(0, pdb_filename.size()-3) + "mtz";

  vector<Cluster> clusters;
  try {
    clusters = find_clusters(pdb_filename, mtz_filename, config);
    if (config.write_1_1_map_to_file) {
      DensityMap density_map;
      density_map.read_from_mtz(mtz_filename, "DELFWT", "PHDELWT");
      density_map.write_ccp4(config.write_1_1_map_to_file);
    }
  }
  catch (const std::runtime_error &e) {
    err(e.what());
  }
  catch (const clipper::Message_fatal &e) {
    err(e.text());
  }

  printf("%ld clusters (with given criteria) found\n", clusters.size());
  for (size_t i = 0; i < clusters.size(); ++i) {
    const clipper::Vec3<>& trn = clusters[i].loc.trn();
    printf("#%-4ld %-3ld grid points, score %-8.4g  (%7.2f,%7.2f,%7.2f)\n",
           i, clusters[i].points.size(), clusters[i].score,
           trn[0], trn[1], trn[2]);
  }

  return 0;
}

#endif // BUILD_EXE

// vim: expandtab:ts=2:sw=2

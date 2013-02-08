
#include <string>
#include <vector>
class CMMDBManager;
namespace clipper { template<class T> class Xmap; }

CMMDBManager* read_pdb(const std::string& pdb_name);

struct Blob
{
    double x;
    double y;
    double z;
    int n_points;
    double score;
};

std::vector<Blob> find_blobs(const std::string& pdb_filename,
                             const std::string& mtz_filename,
                             double sigma_level=1.,
                             double min_volume=11.,
                             double mask_radius=2.);

class DensityMap
{
public:
  DensityMap() : xmap_(NULL) {}
  ~DensityMap();
  void read_from_mtz(const std::string& mtz_name,
                     const std::string& f_col, const std::string& phi_col);
  void write_ccp4(const std::string& fname) const;
  double calculate_stddev() const;
  void mask_by_atoms(CMMDBManager *mol, double mask_radius);
  const clipper::Xmap<float>* xmap() const { return xmap_; }
private:
  clipper::Xmap<float>* xmap_;
};


@setlocal
@set testdir=%~dp0

@echo      === Run dimple on artificial thaumatin example ===
@set PYTHONPATH=%testdir%..\..
@call %CCP4%\..\Python2.7\python "%testdir%..\main.py" "%testdir%thaumatin.mtz" "%testdir%thaumatin.pdb" "%testdir%out" %*
@pause

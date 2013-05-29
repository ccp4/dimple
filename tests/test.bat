@setlocal
@set testdir=%~dp0

@echo "Run dimple on artificial thaumatin example"
@set PYTHONPATH=%testdir%..;%PYTHONPATH%
@ccp4.python %testdir%..\dimple.py %testdir%\thaumatin.mtz %testdir%\thaumatin.pdb %testdir%\out %*
pause

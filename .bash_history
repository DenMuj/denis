showq
qstat -f 637812
showq
python -m venv
module load python/3.6.5 
python -m venv
module purge
module load anaconda/3-5.1.0 
python -m venv
showq
ll
module avail
ls
showq
module purge
showq
module load anaconda/3-5.1.0 
pip list
ls
ll
showq
cd bec3d/
module load git/2.9.5 
git pull
cd bec3d/
cat Par.py 
cd MojOutput/
ll
rm *.err
rm *.out
ll
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
car real3d-rms.txt 
ll
cat 638465.paradox.ipb.ac.rs.err
rm *.err
rm *.out
git add .
git commit -m "malo sutra"
git pull
qsub /home/denis/bec3d/job.pbs 
ll
showq
cd bec3d/MojOutput/
ll
show
showq
cd bec3d/MojOutput/
cd plots/
ll
rm *.png
ll
cd ..
ll
showq
cat /home/denis/bec3d/Par.py 
car /home/denis/bec3d/input/real3d-input 
cat /home/denis/bec3d/input/real3d-input 
mkdir python3.11
ll
cd python3.11
curl https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tar.xz
1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;1;112;112;1;0x1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c1;2c
ll
curl -O https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tar.xz
tar -xf Python-3.11.3.tar.xz
ll
cd Python-3.11.3
./configure --prefix=$(pwd)
./configure --enable-optimizations
make
make install
python --version
cd ..
cd ,,
cd ..
module load Python-3.11.3
module load /home/denis/python3.11/Python-3.11.3
ll
cd python3.11/
ll
cd Python-3.11.3
ll
make install
gcc --version
module avail
module load c++/gnu/9.3.0  
make
cd env
source env/bit/activate
source /env/bit/activate
source /bit/activate
source env/bin/activate
source env/bin/activatce
cd ..
source env/bin/activate
cd env/
curl -O https://files.pythonhosted.org/packages/84/a9/2bf119f3f9cff1f376f924e39cfae18dec92a1514784046d185731301281/scipy-1.10.1.tar.gz
tar -xf scipy-<version>.tar.gz
tar -xf scipy-1.10.1.tar.gz
ll
cd scipy-1.10.1
ll
python3.11 setup.py install
cd ..
ll
rm -r scipy-1.10.1.tar.gz 
curl -o https://files.pythonhosted.org/packages/2c/d4/590ae7df5044465cc9fa2db152ae12468694d62d952b1528ecff328ef7fc/numpy-1.24.3.tar.gz
curl -0 https://files.pythonhosted.org/packages/2c/d4/590ae7df5044465cc9f
ll
module load git/2.9.5 
cd bec3d/MojOutput/
rm *.err
rm *.out
git pull
cd bec3d/MojOutput/
cat /home/denis/bec3d/Par.py 
qsub /home/denis/bec3d/job.pbs 
ll
cd plots/
ll
cd ..
ll
showq
ll
which openssl
rm -r env
ll
cd python3.11/
cd Python-3.11.3/
ll
make clean
./configure --enable-optimizations --with-ensurepip=install --with-openssl=/usr/bin/openssl
make -j 4
make clean
module list
module load c++/gnu/9.3.0 
./configure --enable-optimizations --with-ensurepip=install --with-openssl=/usr/bin/openssl
make -j 4
make clear
make clean
cd ..
ll
mkdir openssl
ll
cd openssl
curl -O https://www.openssl.org/source/openssl-3.1.0.tar.gz
tar -xf openssl-3.1.0.tar.xz
tar -xf openssl-3.1.0.tar.gz
ll
rm -r openssl-3.1.0.tar.gz 
cd openssl-3.1.0/
ll
cat IN
cat INSTALL.md 
./config --prefix=/path/to/installation/directory
./config --prefix=/home/denis/python3.11/openssl/
make
cat INSTALL.md 
ll
cd python3.11/
cd openssl/
ll
cd openssl-3.1.0/
ll
/Configure --prefix=/home/denis/python3.11/openssl/openssl-3.1.0/ --openssldir=/home/denis/python3.11/openssl/openssl-3.1.0/ssl
./Configure --prefix=/home/denis/python3.11/openssl/openssl-3.1.0/ --openssldir=/home/denis/python3.11/openssl/openssl-3.1.0/ssl
module list
module load c++/gnu/9.3.0 
./Configure --prefix=/home/denis/python3.11/openssl/openssl-3.1.0/ --openssldir=/home/denis/python3.11/openssl/openssl-3.1.0/ssl
make
ll
make
make clean
./Configure LIST
cd ..
ll
rm -r openssl/
showq
cd python3.11/
cd ..
python3.11 -m venv evn
cd python3.11/
cd Python-3.11.3/
ll
cd ..
ll
cd python3.11/
cd Python-3.11.3
module load c++/gnu/9.3.0  
make install
python --version
ll
cd ..
ll
python3 -m venv env
ll
source env/bin/activate
ll
cd python3.11/
ll
rm *.xz
ll
mkdir libr
ll
cd libr
pip install scipy
python --version
rm Python-3.11.3
cd ..
rm Python-3.11.3/
rm Python-3.11.3/
rm -r  Python-3.11.3/
ll
python --version
pip install --upgrade python
pip install --upgrade Python3.11.3
pip install --upgrade python
pip install --upgrade python==3.11.3
python --V
python --version
deactivate
ll
cd ..
ll
rm -r env
ll
python --version
cd python3.11/
curl -O https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tar.xz
tar -xf Python-3.11.3.tar.xz
ll
cd Python-3.11.3
./configure --enable-optimizations --with-ensurepip=install
module list
make 
ll
make altinstall
make clean
./configure --enable-optimizations --with-ensurepip=install --prefix=$(pwd)
cat /proc/cpuinfo
make -j 4
make altinstall
ll
cd bin
ll
cd ..
cd
ll
python3.11 -m venv env
python --version
cd python3.11/
ll
rm Python-3.11.3.tar.xz
cd Python-3.11.3/
ls
ls -l python
cd ..
ls
cd Python-3.11.3/
python
ls
cd bin
ls
./python3.11 -m venv /home/denis/env
cd ..
ll
source ./env/bin/activate
cd py
cd python3.11/
cd Python-3.11.3/
ll
cd ..
ll
cd env
ll
pip3.11
pip --version
pip3.11 install scipy
cd /home/denis/env/lib/python3.11/site-packages/
ll
cd ..
ll
cd ..
ll
cd ..
ll
cd ..
ll
cd python3.11/
cd Python-3.11.3/
ll
cd bin
ll
cd ..
/home/denis/python3.11/Python-3.11.3/bin/pip3.11 install scipy
pip
pip --version
python --version
/home/denis/env/lib/python3.11/site-packages/pip install scipy
/home/denis/python3.11/Python-3.11.3/bin/pip3.11 download scipy
--user
deactivate
ll
python3 -m venv env1
sourde env1/bin/activate
ll
source env1/bin/activate
ll
cd env1
ll
cd bin
ll
cd ..
pip3.4 install scipy
cd ..
rm env1
rm -r  env1
rm -r env1
ll
cd /home/denis
ll
rm -r env1
ll
dpkg -s libssl-dev openssl
rpm -q libssl-devel openssl
module load python/3.6.5 
python3.6 -m venv env
ll
source ./env/bin/activate
python --version
pip
pip --version
cd env/
ll
pip install scipy
pip install --upgrade pip
pip list
cd ..
ll
rm -r python3.11
deactivate
module purge
source ./env/bin/activate
python --version
pip list
pip install matplotlib.pyplot
pip install matplotlib
pip list
deactivate
source ./env/bin/activate
cd bec3d/MojOutput/
ll
showq
qdel 638476
qsub /home/denis/bec3d/job.pbs 
ll
rm *.err
rm *.out
ll
source /home/denis/env/bin/activate
ll
cd env/
curl -O https://files.pythonhosted.org/packages/2c/d4/590ae7df5044465cc9fa2db152ae12468694d62d952b1528ecff328ef7fc/numpy-1.24.3.tar.gz
tar -xf numpy-1.24.3.tar.gz
ll
rm -r numpy-1.24.3.tar.gz 
pip list
cd numpy-1.24.3/
ll
cat README.md 
ll
cat INSTALL.rst 
module load c++/gnu/9.3.0 
python setup.py build
cd ..
ll
pip --version
pip install numpy
cd ..
deactivate
python3.6 -m venv env1
python3 -m venv env1
ll
source env1/bin/activate
cd env1
pip --version
pip install numpy
pip intall --upgrade pip
pip install --upgrade pip
pip --version
pip
pip install numpy
deactivate
cd ..
rm -r env1
ll
showq
cd bec3d/MojOutput/
ll
cat 638466.paradox.ipb.ac.rs.out
cd bec3d/MojOutput/
git add .
git commit -m "DRUGODRUGO"
git pull
qsub /home/denis/bec3d/Par.py 
ll
rm Par.py.e638478
rm Par.py.o638478 
ll
module avail
git pull
qsub /home/denis/bec3d/job2.py 
ll
cat 638479.paradox.ipb.ac.rs.err
rm *.err
rm *.out
mkdir drugo
ll
source /home/denis/env/bin/activate
cd drugo/
qsub /home/denis/bec3d/job2.py 
ll
cat 638480.paradox.ipb.ac.rs.err
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
cat 638481.paradox.ipb.ac.rs.err
pip list
git pull
ll
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
cat 638482.paradox.ipb.ac.rs.err
git pull
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
cat 638483.paradox.ipb.ac.rs.err
cat 638483.paradox.ipb.ac.rs.out
rm *.err
rm *.out
cd ..
ll
cd drugo/
git pull
qsub /home/denis/bec3d/job2.py 
ll
cat 638484.paradox.ipb.ac.rs.err
ll
rm *.err
rm *.out
rm *.txt
cd ..
rm drugo
rm -r drugo
ll
which openssl
which openssl-devel
yum list installed openssl-devel
rpm -q --queryformat '%{INSTALLPREFIX}\n' openssl-devel
rpm -qi openssl-devel
ls /usr
ls /usr/include/openssl
ls /usr/lib/
which openssl-devel
showq
cd bec3d/MojOutput/
ll
cat 638477.paradox.ipb.ac.rs.err
rm *.err
rm *.out
cd bec3d/MojOutput/
qsub /home/denis/bec3d/job2.py 
ll
cat reald3d-rms.txt
cat real3d-rms.txt 
ll
showq
qdel 638485
module load git/2.9.5 
git add .
git commit -m "LASDLASD"
git pull
cd bec3d/
git add .
git commit -m "aksfaskmfmaskfas"
git pull
module load git/2.9.5 
cd bec3d/
cat Par.py
cat par2.py
cd MojOutput/
ll
rm *.err
rm *.out
qsub /home/denis/bec3d/job2.py 
ll
module list
ll
cat real3d-rms.txt 
ll
cat 638486.paradox.ipb.ac.rs.out
git pull
git add .
git commit -m "KAKKSDA"
git pull
cd bec3d/input/real3d-input 
cat bec3d/input/real3d-input 
git add .
git commit -m "WWEQWEQ"
git pull
module load git/2.9.5 
git pull
git add .
git commit -m "AsDASDw"
git pull
cd bec3d/input/real3d-input 
cd bec3d/
git add .
git commit -m "ASDASDWwq"
git pull
cd bec3d/
cd MojOutput/
ll
rm *.err
rm *.out
git pull
git add .
git commit -m "WCHANGES"
git pull
cat /home/denis/bec3d/input/real3d-input 
git pull
module load git/2.9.5 
git pull
cat /home/denis/bec3d/input/real3d-input 
nano /home/denis/bec3d/input/real3d-input 
cat /home/denis/bec3d/input/real3d-input 
nano /home/denis/bec3d/par2.py 
cat /home/denis/bec3d/par2.py 
nano /home/denis/bec3d/Par.py 
qsub /home/denis/bec3d/job2.py 
ll
cat real3d-rms.txt 
ll
cat real3d-rms.txt 
cd bec3d/MojOutput/
ll
showq
qdel 638487
rm *.err
rm *.out
cat /home/denis/bec3d/Par.py
qsub /home/denis/bec3d/job.pbs 
ll
showq
cat real3d-rms.txt 
cat /home/denis/bec3d/Par.py 
showq
ll
source ./env/bin/activate
python --version
import ssl
which openssl
import ssl
openssl --version
showq
cd bec3d/
module load git/2.9.5 
git add .
git commit -m "OPSCIP"
git pull
cd bec3d/MojOutput/
git pull
git add .
git commit -m "AWWWQQ"
git pull
cd ..
cat par2.py
git pull
cd MojOutput/
qsub /home/denis/bec3d/job2.py 
ll
rm 638488.paradox.ipb.ac.rs.out
rm 638488.paradox.ipb.ac.rs.err
ll
showq
cd bec3d/MojOutput/
ll
cat 638488.paradox.ipb.ac.rs.out
cd bec3d/MojOutput/
ll
cd bec3d/MojOutput/
ll
cat real3d-rms.txt 
showq
cat /home/denis/bec3d/input/real3d-input 
cat /home/denis/bec3d/MojOutput/real3d-rms.txt 
cat /home/denis/bec3d/Par.py 
cd /home/denis/bec3d/MojOutput/
qsub /home/denis/bec3d/job.pbs 
qdel 638575
ll
rm 638575.paradox.ipb.ac.rs.err
rm 638575.paradox.ipb.ac.rs.out
git push
git pull
qsub /home/denis/bec3d/job.pbs 
ll
showq
cd /home/denis/bec3d/MOj
cd /home/denis/bec3d/MojOutput/
ll
cat 638523.paradox.ipb.ac.rs.out
showq
cd /home/denis/bec3d/MojOutput/
ll
showq
cat /home/denis/bec3d/MojOutput/real3d-rms.txt 
showq
cd /home/denis/bec3d/MojOutput/
ll
cat 638576.paradox.ipb.ac.rs.out
rm 638576.paradox.ipb.ac.rs.out
rm 638576.paradox.ipb.ac.rs.err
cd bec3d/
module load git/2.9.5 
git pull
cd bec3d/MojOutput/
module load git/2.9.5 
git pull
cd bec3d
cat Par.py 
cd MojOutput/
qsub /home/denis/bec3d/job.pbs 
ll
showq
qstat -f 637996
cd bec3d/MojOutput/
ll
car /home/denis/bec3d/input/real3d-input 
cat /home/denis/bec3d/input/real3d-input 
showq
qstat -f 638932
cd bec3d/MojOutput/
module load git/2.9.5 
git add .
git commit -m "SA AMPLTIR"
git pull
cd bec3d/
git pull
git add .
git commit -a
cd bec3d/MojOutput/
module load git/2.9.5 
git add .
git commit -m "AS@2r2A"
git push
git pull
cd bec3d/MojOutput/
git status
git commit -m "ADSASwtw"
git status
git push origin master
git status
git pill
git pull
cat /home/denis/bec3d/Par.py 
qsub /home/denis/bec3d/job.pbs 
ll
kk
ll
showq
ll
cat 638985.paradox.ipb.ac.rs.err
rm *.err
rm *.out
qsub /home/denis/bec3d/job.pbs 
ll
;;
ll
cd bec3d/MojOutput/
ll
showq
module load git/2.9.5 
git pull
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
rm *.err
rm *.out
git pull
qsub /home/denis/bec3d/job.pbs 
ll
cat /home/denis/bec3d/Par.py 
ll
showq
ll
showq
ll
showq
ll
showq
qdel 638988
ll
rm *.err
rm *.out
git pull
qsub /home/denis/bec3d/job.pbs 
ll
cd bec3d/MojOutput/
ll
cat 638848.paradox.ipb.ac.rs.out
ce bec3d/MojOutput/
cd bec3d/MojOutput/
showq
vi
cd bec3d/MojOutput/o
cd bec3d/MojOutput/
ll
cat638989.paradox.ipb.ac.rs.err
cat 638989.paradox.ipb.ac.rs.err
vi /home/denis/bec3d/Par.py 
cat /home/denis/bec3d/Par.py 
vi /home/denis/bec3d/Par.py 
cat /home/denis/bec3d/Par.py 
cd bec3d/
module load git/2.9.5 
git add .
git commit -m "PARAGON"
git status
git push
git status
cd MojOutput/
ll
rm *.err
rm *.put
rm *.out
qsub /home/denis/bec3d/job.pbs 
ll
l
ll
nvim
ll
nvi
vi
vim
ll
showq
cd bec3d/MojOutput/
ll
showq
cd bec3d/MojOutput/
ll
module load git/2.9.5
git pull
ll
cat /home/denis/bec3d/Par.py 
qsub /home/denis/bec3d/job.pbs 
ll
639011.paradox.ipb.ac.rs.out
cat 639011.paradox.ipb.ac.rs.out
ll
rm 639011.paradox.ipb.ac.rs.err
showq
ll
cd bec3d/MojOutput/
ll
cat real3d-rms.txt 
ll
cd bec3d/MojOutput/
xeyes
cd bec3d/MojOutput/
xeyes
cd bec3d/MojOutput/
xeyes
cd plots/
ll
display plot_40.png 
cd bec3d/MojOutput/plots/
ll
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/
cd bec3d/MojOutput/plots/
cd plots/
display plot_40.png 
echo $DISPLAY
cd ..
ll
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_40.png 
cd bec3d/MojOutput/plots/
ll
display plot_41.png 
cd bec3d/
module load git/2.9.5 
git status
showq
cd MojOutput/plots/
ll
cd ..
git add .
git commit -m "PLOTS"
git push
ll
cd MojOutput/plots/
cd bec3d/MojOutput/plots/
ll
display plot_42.png 
display plot_42.png,plot_41.plot 
display plot_41.png 
display plot_40.png 
display plot_42.png 
cd bec3d/MojOutput/plots/
ll
display plot_40.png 
cd bec3d/MojOutput/plots/
display plot_43.png 
display plot_44.png 
display plot_45.png 
display plot_46.png 
display plot_47.png 
display plot_48.png 
display plot_49.png 
display plot_50.png 
display plot_51.png 
:ll
ll
cd bec3d/MojOutput/plots/
ll
showq
top
cd bec3d/MojOutput/
ll
cd bec3d/MojOutput/
ll
cat real3d-rms.txt 
cd bec3d/MojOutput/
ll
cat real3d-rms.txt 
cd plots/
ll
cd ..
ll
cat real3d-rms.txt 
ll
cat 639019.paradox.ipb.ac.rs.out
cd bec3d/MojOutput/
module load git/2.9.5 
git add .
git commit -m "SLIKEEE"
git push
showq
cd bec3d/
module load git/2.9.5 
git status
git pull
cd MojOutput/
qsub /home/denis/bec3d/job.pbs 
ll
showq
ll
cd bec3d/MojOutput/
ll
cat 639370.paradox.ipb.ac.rs.out
cd plots/
ll
cd ..
git add .
git commit -m "ASLFWWW"
git push
cd MojOutput/
ll
cat 639370.paradox.ipb.ac.rs.out
showq
qstat -f 639601
showq
qstat -f 645158
qstat -f 645063
showq
ls -la
joe .bash_profile 
exit
module list
ll
tar -xvzf denis.tar.gz 
ll
cd denis
ll
cat makefile 
make imre3d-ms-ddiX-rot-mpi-qf-grad-mu-self-Q35s
ll
rm *.o
ll
cd ..
ll
mv denis DBEC
ll
mkdir test
cd test/
joe imag3d.pbs
ll
joe input 
rm *~
ll
showq
qsub imag3d.pbs 
ll
cat imag3d-mu.txt 
cat imag3d-out.txt 
cat imag3d-rms.txt 
cat imag3d-mu.txt 
cat imag3d-rms.txt 
joe ../.bash_profile
. ../.bash_profile
bin2vtk_yz-aho 
bin2vtk_yz-aho imag3d-den-niter 10000 1 256 256 256 1
cat imag3d-mu.txt 
bin2vtk_yz-aho imag3d-den-niter 10000 1 256 256 256 1
bin2vtk_xy-aho imag3d-den-niter 10000 1 256 256 256 1
cat imag3d-rms.txt 
showq
qdel 649158
ll
ll 6*
rm 6* *bin *txt *vtk
ll
joe input 
rm *~
qsub imag3d.pbs 
cat imag3d-out.txt 
cat imag3d-rms.txt 
showq
qstat -f 649167
ssh gn021
su -
cat imag3d-rms.txt 
bin2vtk_xy-aho imag3d-den-niter 10000 1 256 256 256 1
bin2vtk_yz-aho imag3d-den-niter 10000 1 256 256 256 1
bin2vtk_xy-aho imag3d-den-niter 10000 1 256 256 256 1
showq
qdel 649167
ll
showq
ll 6*
showq
ls -la
rm -r denis.tar.gz 
ls -la
cd test
ll
ls -lh
git status
module avail
module load git/2.9.5 
cd test
git add .
git push
cd ..
git add .
git push
git push test
ll
git status
git push

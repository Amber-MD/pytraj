pacman -Sy pacman 
pacman -Syu
pacman -Su
pacman -S unzip make rsync git tcsh libbz2-devel zilb-devel diffutils tar
pacman -S --noconfirm --needed mingw-w64-x86_64-openblas
pacman -S --noconfirm --needed  mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-ncurses mingw-w64-x86_64-readline

# python
pacman -S python2
pacman -S mingw-w64-x86_64-python2-numpy
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py

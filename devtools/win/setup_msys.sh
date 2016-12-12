# pacman -Sy --noconfirm --needed pacman 
# pacman -Syu --noconfirm --needed
# pacman -Su --noconfirm --needed
pacman -S --noconfirm --needed unzip make rsync git tcsh libbz2-devel zilb-devel diffutils tar
pacman -S --noconfirm --needed mingw-w64-x86_64-openblas
pacman -S --noconfirm --needed  mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-ncurses mingw-w64-x86_64-readline

# python
pacman -S --noconfirm --needed python2
pacman -S --noconfirm --needed mingw-w64-x86_64-python2-numpy
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py

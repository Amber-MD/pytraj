
if [ ! -f /usr/local/lib/libnetcdf.a ]; then
  curl -fsS -o netcdf-4.3.3.zip ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.zip
  # unzip netcdf-4.3.3.zip
  7z x netcdf-4.3.3.zip
  cd netcdf-4.3.3
  exec 0</dev/null
  ./configure --enable-static --disable-netcdf-4 --prefix=/usr/local/ --disable-dap
  make -r install;
else
  echo 'Have Cached NetCDF';
fi

appveyor DownloadFile https://ci.appveyor.com/api/buildjobs/hfs561h7ex0r4mtp/artifacts/cpptraj-6082574.zip -FileName cpptraj.zip
mkdir cpptraj
mv cpptraj.zip cpptraj/
cd cpptraj/
7z x cpptraj.zip
ls lib/
cd ../

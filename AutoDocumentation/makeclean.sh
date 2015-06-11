mv source/index.rst ./index.rst
rm -rf source/*.rst
mv ./index.rst source/index.rst
make clean
sphinx-apidoc -o ./source ../
make html

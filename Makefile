PYTHON="/home/bkille/packages/debug_python/bin/python3.7"


bitarray/_bitarray.so: bitarray/_bitarray.cpp
	$(PYTHON) setup.py build_ext --inplace


test: bitarray/_bitarray.so
	$(PYTHON) -c "import bitarray; bitarray.test()"


doc: bitarray/_bitarray.so
	$(PYTHON) update_readme.py
	$(PYTHON) setup.py sdist
	twine check dist/*


clean:
	rm -rf build dist
	rm -f bitarray/*.o bitarray/*.so
	rm -f bitarray/*.pyc
	rm -f examples/*.pyc
	rm -rf bitarray/__pycache__ *.egg-info
	rm -rf examples/__pycache__

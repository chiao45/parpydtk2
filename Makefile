include make.inc

.pyx.cpp:
	$(CYTHON) --cplus -$(PYTHON_VERSION) -I./pydtk2 -X boundscheck=False $<

inplace:
	env CC=$(CXX) $(PYTHON) setup.py build_ext --inplace

.PHONY: inplace

install:
	env CC=$(CXX) $(PYTHON) setup.py install

.PHONY: install

regen_cython: $(PYX_OBJ)

clean_cython:
	rm -rf $(PYX_OBJ)

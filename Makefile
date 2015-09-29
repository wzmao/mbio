# Makefile for some developer commands

REPOPATH = `pwd`

.PHONY: help build build3 remove test

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  build		to build extensions in place"
	@echo "  remove 	to remove contributed modules"

build2.6:
	python2.6 setup.py build_ext --inplace --force

build2:
	python2 setup.py build_ext --inplace --force

build:
	python setup.py build_ext --inplace --force

build3:
	python3 setup.py build_ext --inplace --force

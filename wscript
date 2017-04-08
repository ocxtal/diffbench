#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.recurse('gaba')
	opt.load('compiler_c')
	opt.load('compiler_cxx')

def configure(conf):
	conf.options.bit = 2
	conf.recurse('gaba')

	conf.load('ar')
	conf.load('compiler_c')
	conf.load('compiler_cxx')

	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')

	conf.env.append_value('CXXFLAGS', '-O3')
	conf.env.append_value('CXXFLAGS', '-std=c++11')
	conf.env.append_value('CXXFLAGS', '-march=native')

	conf.env.append_value('LIB_DIFFBENCH', conf.env.LIB_GABA)
	conf.env.append_value('OBJ_DIFFBENCH', ['aed.o', 'alinear.o', 'aaffine.o', 'edlib.o', 'bench.o'] + conf.env.OBJ_GABA)

def build(bld):
	bld.recurse('gaba')

	bld.objects(source = 'aed.c', target = 'aed.o')
	bld.objects(source = 'alinear.c', target = 'alinear.o')
	bld.objects(source = 'aaffine.c', target = 'aaffine.o')
	bld.objects(source = 'edlib.cpp', target = 'edlib.o')
	bld.objects(source = 'bench.c', target = 'bench.o', defines = ['BENCH'])

	bld.program(source = 'empty.cpp', target = 'bench', use = bld.env.OBJ_DIFFBENCH, lib = bld.env.LIB_DIFFBENCH)

	"""
	bld.program(
		target = 'bench',
		use = bld.env.OBJ_DIFFBENCH,
		lib = bld.env.LIB_DIFFBENCH,
		defines = ['BENCH'])
	"""

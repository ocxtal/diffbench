#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.recurse('gaba')
	opt.load('compiler_c')

def configure(conf):
	conf.options.bit = 2
	conf.recurse('gaba')

	conf.load('ar')
	conf.load('compiler_c')

	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')

	conf.env.append_value('LIB_DIFFBENCH', conf.env.LIB_GABA)
	conf.env.append_value('OBJ_DIFFBENCH', ['aed.o', 'alinear.o', 'aaffine.o'] + conf.env.OBJ_GABA)

def build(bld):
	bld.recurse('gaba')

	bld.objects(source = 'aed.c', target = 'aed.o')
	bld.objects(source = 'alinear.c', target = 'alinear.o')
	bld.objects(source = 'aaffine.c', target = 'aaffine.o')

	bld.program(
		source = ['bench.c'],
		target = 'bench',
		use = bld.env.OBJ_DIFFBENCH,
		lib = bld.env.LIB_DIFFBENCH,
		defines = ['BENCH'])

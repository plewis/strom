project('strom', 'cpp',
	default_options : ['cpp_std=c++11'],
	version : '3.0-dev')
dep_beagle = dependency('hmsbeagle-1', version : '>=3')
# Boost doesn't ship with pkg-config *.pc files per ticket 1094
# https://svn.boost.org/trac10/ticket/1094, so we have to explicitly import the
# required library using meson's compiler.find_library().
compiler = meson.get_compiler('cpp')
dep_boost = compiler.find_library('boost_program_options')
dep_eigen = dependency('eigen3')
executable('strom', 'main.cpp',
	   dependencies : [dep_beagle,
			   dep_boost,
			   dep_eigen])

project('strom', 'cpp',
	default_options : ['cpp_std=c++11'],
	version : '3.0-dev')
dep_beagle = dependency('hmsbeagle-1', version : '>=3')
dep_boost = dependency('boost',
                       modules : ['program_options'])
dep_ncl = dependency('nclv2.1')
dep_eigen = dependency('eigen3')
executable('strom', 'main.cpp',
	   dependencies : [dep_beagle,
			   dep_boost,
			   dep_ncl,
			   dep_eigen])

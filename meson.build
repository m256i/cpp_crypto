project('cpp_crypto', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20'])

exe = executable(
  'cpp_crypto', 
  ['main.cpp', 'bigint.cpp'],
  link_args : ['-lgmp'], 
  cpp_args: ['-g', '-O3', '-mtune=native', '-march=native'], 
  install : true)

test('basic', exe)

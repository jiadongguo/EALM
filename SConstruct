import os,sys
# 定义程序生成的各个输出文件的名字
cstd_lib='cstd' #自定义库文件的名字
main_model='jd_model'
main_ricker='jd_ricker'
# 自身库文件
comm_libs='''
m openblas
'''
env = Environment(
    CC='gcc',  # 指定 C 编译器
    CCFLAGS=['-w','-fPIC', '-g','-O3'], 
    SHLIBFLAGS=['-shared'],
    LIBPATH=['.']
)

# 生成标准库文件
cstd_src='''
cstd.c sim.c
'''
env.SharedLibrary(
    target=cstd_lib, 
    source=Split(cstd_src),
    LIBS=Split(comm_libs)
)
# 雷克子波生成
env.Program(
    target=main_ricker, 
    source=['ricker.c'],
    LIBS=Split(cstd_lib+comm_libs)
)
# main_model的源文件和库文件
model_src='''
main_model.c a2d_mod_28.c abc2d.c
'''

program = env.Program(
    target=main_model, 
    source=Split(model_src),
    LIBS=Split(cstd_lib+comm_libs)
)
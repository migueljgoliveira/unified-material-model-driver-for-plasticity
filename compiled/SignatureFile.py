
def main():

    s = " "

    file = open('ummdp_vfm.pyf','r')
    pyf = file.readlines()
    file.close()

    npyf = pyf.copy()
    i = 0
    for line in pyf:
        try:
            if line.strip(' ').split()[0] == 'subroutine':
                npyf.insert(i+1,f'{s*12}use f2py_stop__user__routines\n')
                npyf.insert(i+2,f'{s*12}intent(callback,hide) :: f2py_stop\n')
                npyf.insert(i+3,f'{s*12}external f2py_stop\n')

                i += 3
        except:
            pass

        i += 1

    npyf.insert(0+3,'python module f2py_stop__user__routines\n')
    npyf.insert(1+3,f'{s*4}interface f2py_stop_user_interface\n')
    npyf.insert(2+3,f'{s*8}subroutine f2py_stop()\n')
    npyf.insert(3+3,f'{s*12}intent(callback,hide) f2py_stop\n')
    npyf.insert(4+3,f'{s*8}end subroutine f2py_stop\n')
    npyf.insert(5+3,f'{s*4}end interface f2py_stop_user_interface\n')
    npyf.insert(6+3,'end python module f2py_stop__user__routines\n')
    npyf.insert(7+3,'\n')

    pyf = open('ummdp_vfm2.pyf','w')
    pyf.writelines(npyf)

if __name__ == '__main__':

    main()
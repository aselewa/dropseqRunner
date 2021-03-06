
#!/usr/bin/env python

import filecmp
import config

def compare_dir(dir1, dir2, ignore=None):

    dir_cmp = filecmp.dircmp(dir1, dir2, ignore=ignore)
    if dir_cmp.left_only:
        print('LEFT ONLY FILES: ', dir_cmp.left_only)
        return False
    if dir_cmp.right_only:
        print('RIGHT ONLY FILES: ', dir_cmp.right_only)
        return False

    _, mismatch, errors = filecmp.cmpfiles(dir1, dir2, dir_cmp.common_files, shallow=False)

    if mismatch:
        print('MISMATCH: ', mismatch)
        return False
    if errors:
        print('ERRORS: ', errors)
        return False
    
    return True

if __name__ == "__main__":
    print("\n\n")
    print("**********************************************")
    assert compare_dir(config.expected_dir(), config.mtx_dir(), ignore=['.DS_Store']), 'outputs are NOT the same!. Test workflow did not run properly. Make sure your conda environment was created using the given environment file.\n*********************************************\n\n'
    print('Yay! Outputs are identical and test workflow ran properly.')
    print("**********************************************\n\n")

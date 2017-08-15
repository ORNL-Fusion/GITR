#python library tools for gitr

from distutils.dir_util import copy_tree

def copy_folder(from_folder, to_folder):
    copy_tree(from_folder, to_folder)


if __name__ == "__main__":
    print 'no functions'

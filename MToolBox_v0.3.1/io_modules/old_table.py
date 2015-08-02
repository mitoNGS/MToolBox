# encoding=utf8

# import sys
# if '..' not in sys.path:
#     sys.path.append('..')

from classifier import tree

def write_haplogroup(file_handle, rm_file, haplogroup):
    #solo quelli filtrati
    filtered_positions = tree.filter_positions(haplogroup)
    print "CIAO"
    print haplogroup
    print filtered_positions
    pos_list = []
    for position in filtered_positions:
        pos_list.append( (haplogroup.name,) + position.print_table() )
    file_handle.writerows(pos_list)

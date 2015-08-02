# coding=utf8
"""
Legge i file contenuti nella directory e inizializza i formati di file
supportati in lettura/scrittura
"""

import os, imp

plugins = set()
for plugin in os.listdir(os.path.join(os.path.dirname(__file__),'.')):
    if plugin.find('.py') != -1:
        plugins.add(plugin[:plugin.find('.py')])
plugins.remove('__init__')

write_funcs = {}
"""
Variabile contenente le funzioni per scrivere il formato
"""
load_funcs = {}
"""
Variabile contenente le funzioni per leggere il formato
"""
check_funcs = {}
"""
Variabile contenente le funzioni per individuare il formato
"""
file_types = {}
"""
Variabile contenente il tipo di formato (0, liscio, 1, complesso)
"""
func_types = {'write': write_funcs, 'load': load_funcs, 'check': check_funcs}
"""
Variabile contenente i puntatori ai dizionari
"""

for plugin in plugins:
    fp, pathname, description = imp.find_module(plugin, __path__)#['.','..', 'files'])
    try:
        mod = imp.load_module(plugin, fp, pathname, description)
    finally:
        if fp:
            fp.close()
    file_types[plugin] = getattr(mod, 'FILETYPE')
    for func_type in func_types:
        try:
            func_types[func_type][plugin] = getattr(mod, func_type+'_'+plugin)
        except:
            pass

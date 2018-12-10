#!usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#referencia caso precise na documentação
#https://www.python.org/doc/essays/graphs/
#Copyright (c) 1998, 2000, 2003 Python Software Foundation.
#All rights reserved.
#Licensed under the PSF license.
#
#
#
class GraphDAO (object):

    #se funcionar a solução á daqui:https://www.geeksforgeeks.org/generate-graph-using-dictionary-python/
    def create_graph(self,lista_de_alelos,ultima_posicao):
        grafo = {}
        tupla = []
        cont=1
        #nao ta mais dando problema
        last_pos = ultima_posicao
        #print (last_pos)
        for i in range(len(lista_de_alelos)):
            #grafo = [lista_de_alelos[i]]
            for j in range(len(lista_de_alelos)):
                if lista_de_alelos[j].snp_pos == cont + 1:
                    tupla.append(lista_de_alelos[j])

            grafo[lista_de_alelos[i]] = tupla
            tupla =[]
            if lista_de_alelos[i].snp_pos != last_pos and lista_de_alelos[i+1].snp_pos == (cont + 1):
                cont = cont + 1
        return grafo

    #método de retornar todos os caminhos
    #método da documentação do python
    def find_all_paths(self,graph, start, end, path=[]):
            path = path + [start.nome]
            if start == end:
                return [path]
            if not graph.get(start):
                return []
            paths = []
            for node in graph[start]:
                if node not in path:
                    newpaths = self.find_all_paths(graph, node, end, path)
                    for newpath in newpaths:
                        paths.append(newpath)
            return paths

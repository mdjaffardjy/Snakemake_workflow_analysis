#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 08:55:23 2021

@author: marinedjaffardjy
"""

import numpy as np
import pandas as pd

# here is a collection of functions that help 

def get_parsing_stats(df_800, dict_tool, total_extr, tot_tools):
    # prints the stats of the parser
    print("total wf parsed : "+ str(len(df_800)))
    res = 0
    for i in range(0, len(df_800["tool names"])):
        if df_800["tool names"][i] != set():
            res += 1
    res
    print("wf with tools : " + str(res))
    print("nb of tools found : " + str(len(total_extr)))
    print("nb of tools found with a corresponding biotools entry : " + str(len(dict_tool)))
    print("nb of individual tools found : " + str(len(tot_tools)))
    print("mean number of tools per wf (only wf with at least one tool) : "+ str(len(dict_tool)/res))


def get_bag_of_toolnames(tools):
    names = []
    for tool in tools:
        names.append(tool['name'])
    names = set(names)
    return names

def get_bag_of_topics(tools):
    #a function that gets all the topic information of the tools in the workflow
    #input : a list of dict of the tools informations (with the toolname, its biotools ID, colletion ID, {topic}, {function})
    #output : a set of the topics names and a set of the topics uris
    topics_uris = []
    topics_labels = []
    #there are several words in 'topics', for now we take 
    #all of them but we could choose
    for tool in tools:
        for element in tool['topic'][0]:
            topics_labels.append(element['term'])
            topics_uris.append(element['uri'])           
    return set(topics_labels), set(topics_uris)

def get_bag_of_topics_syn(tools):
    #a function that gets all the topic information of the tools in the workflow
    #input : a list of dict of the tools informations (with the toolname, its biotools ID, colletion ID, {topic}, {function})
    #output : a set of the topics names and a set of the topics uris
    topics_uris = []
    topics_labels = []
    #there are several words in 'topics', for now we take
    #all of them but we could choose
    for tool in tools:
        for element in tool['topic'][1]:
            topics_labels.append(element['term'])
            topics_uris.append(element['uri'])
    return set(topics_labels), set(topics_uris)

def get_bag_of_operations(tools):
    #a function that gets alll the operations information of the tools in the workflow
    #input : a list of dict of the tools informations (with the toolname, its biotools ID, colletion ID, {topic}, {function})
    #output : a set of the operation names and a set of the operations uris
    op_uris = []
    op_labels = []
    #there are several functions in operations, we consider all of them
    for tool in tools:
        for f in tool['function']:
            for element in f['operation'][0]:
                op_labels.append(element['term'])
                op_uris.append(element['uri'])           
    return set(op_labels), set(op_uris)

def get_bag_of_operations_syn(tools):
    #a function that gets alll the operations information of the tools in the workflow
    #input : a list of dict of the tools informations (with the toolname, its biotools ID, colletion ID, {topic}, {function})
    #output : a set of the operation names and a set of the operations uris
    op_uris = []
    op_labels = []
    #there are several functions in operations, we consider all of them
    for tool in tools:
        for f in tool['function']:
            for element in f['operation'][1]:
                op_labels.append(element['term'])
                op_uris.append(element['uri'])
    return set(op_labels), set(op_uris)

def get_bag_of_inputs(tools):
    #a function that gets all the inputs of the tools in the workflow
    #input : a list of dict of the tools informations (with the toolname, its biotools ID, colletion ID, {topic}, {function})
    #output : a set of the input names and a set of the inputs uris (we consider the data part)
    in_uris = []
    in_labels = []
    #there are several functions in operations, we consider all of them
    if tools != []:
        for tool in tools:
            for f in tool['function']:
                for element in f['input']:

                    if ('data' in element.keys()) :
                        in_labels.append(element['data']['term'])
                        in_uris.append(element['data']['uri'])
                    else :
                        in_labels.append(element['term'])
                        in_uris.append(element['uri'])



            
    return set(in_labels), set(in_uris)

def get_bag_of_outputs(tools):
    #a function that gets all the outputs of the tools in the workflow
    #input : a list of dict of the tools informations (with the toolname, its biotools ID, colletion ID, {topic}, {function})
    #output : a set of the output names and a set of the outputs uris (we consider the data part)
    out_uris = []
    out_labels = []
    #there are several functions in operations, we consider all of them
    in_uris = []
    in_labels = []

    # there are several functions in operations, we consider all of them
    if tools != []:
        for tool in tools:
            for f in tool['function']:
                for element in f['output']:

                    if ('data' in element.keys()):
                        out_labels.append(element['data']['term'])
                        out_uris.append(element['data']['uri'])
                    else:
                        out_labels.append(element['term'])
                        out_uris.append(element['uri'])
            
    return set(out_labels), set(out_uris)

#TODO : faire un seul df/outil avec le nb de wf etc par outil. prendre le nom de l'outil dans le dict et pas dans la liste tot tools
def get_list_total_tools(df):
    total_tools = []
    for tools in df['tool names']:
        total_tools.extend(tools)
    return set(total_tools)

def get_df_attributes(wf_list):
    df = pd.DataFrame(columns=["name","nb of rules", "nb of tools", "tool names","topics", "operations","inputs", "outputs"])

    for wf in wf_list :
        print(wf['filename'])
        #TODO : tqdm : barre de progression
        nb_rules = wf["rulecount"]
        tool_names = get_bag_of_toolnames(wf['tools annotations'])
        nb_tools = len(tool_names)
        topics, x = get_bag_of_topics(wf['tools annotations'])
        topics_syn, top_uris = get_bag_of_topics_syn(wf['tools annotations'])
        operations, op_uris = get_bag_of_operations(wf['tools annotations'])
        operations_syn, x = get_bag_of_operations_syn(wf['tools annotations'])
        inp, x = get_bag_of_inputs(wf['tools annotations'])
        out, x = get_bag_of_outputs(wf['tools annotations'])
        df = df.append({"name": wf['filename'],
                        "nb of rules": nb_rules,
                        "nb of tools": nb_tools,
                        "tool names": tool_names,
                        "topics": topics,
                        "top_uris": top_uris,
                        "topics_syn": topics_syn,
                        "operations": operations,
                        "op_uri" : op_uris,
                        "operations_syn": operations_syn,
                        "inputs": inp,
                        "outputs": out}, ignore_index=True)
    return df

def get_co_occurence(wf,tool1,tool2):
    #inputs : list of tools for a wf, string of tool1 name, str of tool2 name
    #output : the co-occurence matrix of all the tools
    #print(wf)
    #print("1 " + tool1+" 2 "+tool2 )
    occ = 0
    if (tool1 in wf and tool2 in wf):
        occ = 1
    return occ

# TODO : n'en faire qu'un objet ?
def get_df_stats_annot(dict_tools, tot_tools):
    #ne prend pas en compte les synonymes
    df_tools_stats_annot = pd.DataFrame(columns=["toolname",
                                                 "tool_id",
                                                 "nb_operations",
                                                 "nb_inputs",
                                                 "nb_outputs",
                                                 "nb_topics"])

    for tool in dict_tools: #TODO : unicité ? dict tools plus gd que tot tools
        # print(tool)
        if len(dict_tools[tool]["function"]) > 0:
            nb_operations = len(dict_tools[tool]["function"][0]["operation"][0])
            nb_inputs = len(dict_tools[tool]["function"][0]["input"])
            nb_outputs = len(dict_tools[tool]["function"][0]["output"])
        else:
            nb_operations = 0
            nb_inputs = 0
            nb_outputs = 0
        nb_topics = len(dict_tools[tool]["topic"][0])
        df_tools_stats_annot = df_tools_stats_annot.append({"toolname": tool,
                                                            "tool_id" : dict_tools[tool]["name"],
                                                            "nb_operations": nb_operations,
                                                            "nb_inputs": nb_inputs,
                                                            "nb_outputs": nb_outputs,
                                                            "nb_topics": nb_topics}, ignore_index=True)
    return df_tools_stats_annot

def get_df_stats_annot_with_syn(dict_tools, tot_tools):
    #ne prend pas en compte les synonymes
    df_tools_stats_annot = pd.DataFrame(columns=["toolname",
                                                 "tool_id",
                                                 "nb_operations",
                                                 "nb_inputs",
                                                 "nb_outputs",
                                                 "nb_topics"])

    for tool in dict_tools: #TODO : unicité ? dict tools plus gd que tot tools
        # print(tool)
        if len(dict_tools[tool]["function"]) > 0:
            nb_operations = len(dict_tools[tool]["function"][0]["operation"][0])+len(dict_tools[tool]["function"][0]["operation"][1])
            nb_inputs = len(dict_tools[tool]["function"][0]["input"])
            nb_outputs = len(dict_tools[tool]["function"][0]["output"])
        else:
            nb_operations = 0
            nb_inputs = 0
            nb_outputs = 0
        nb_topics = len(dict_tools[tool]["topic"][0])+len(dict_tools[tool]["topic"][1])
        df_tools_stats_annot = df_tools_stats_annot.append({"toolname": tool,
                                                            "tool_id" : dict_tools[tool]["name"],
                                                            "nb_operations": nb_operations,
                                                            "nb_inputs": nb_inputs,
                                                            "nb_outputs": nb_outputs,
                                                            "nb_topics": nb_topics}, ignore_index=True)
    return df_tools_stats_annot


#"TODO : FAIRE UNE FCT POUR ÇA"

def get_df_tools(df, tot_tools):
    df_tools = pd.DataFrame(columns=["toolname","wf","wf_nb"])

    for tool in tot_tools:
        wk=[]
        #print("tool "+tool)
        for i in range(0,len(df["tool names"])):
            if(str(tool) in df["tool names"][i]):
                #print(df["tool names"][i])
                wk.append([i, df["name"][i]])
        #print(str(wk))
        df_tools=df_tools.append({"toolname" : str(tool),
                          "wf" : wk,
                          "wf_nb" : len(wk)}, ignore_index = True)
    return df_tools

#faire un dataframe avec comme labels de lignes et de col les noms d'outils
def get_df_cooc_tools(df_800, tot_tools):
    df_cooc_tools = pd.DataFrame(columns=tot_tools,index = tot_tools)
    for tool1 in tot_tools :
        for tool2 in tot_tools :
            df_cooc_tools[tool1][tool2]=0
            for i in range(0,len(df_800["tool names"])) :
                df_cooc_tools[tool1][tool2] += get_co_occurence(df_800["tool names"][i],tool1,tool2) #soit 0, n si trouvé n fois,
    return df_cooc_tools

def get_proximity_score(wf1,wf2):
    # inputs : list of topics/operations/inputs/outputs for a wf1 and for a wf2
    # output : float of the proximity score
    s1 = set(wf1)
    s2 = set(wf2)
    if s1 != set() and s2 != set():
        return float(len(s1.intersection(s2)) / len(s1.union(s2)))
    else:
        return 0.0

#score de proximité - on peut utiliser les topics et les opérations
def get_prox_dataframe(df, sujet):
    # dataframe de workflows, sujet = string pour le "thème" sur lequel on fait la matrice de proximité
    df_prox = pd.DataFrame(columns=df["name"], index=df["name"])
    for i in range(0, len(df["tool names"])):
        for j in range(0, len(df["tool names"])):
            df_prox[df["name"][i]][df["name"][j]] = get_proximity_score(df[sujet][i].extend(df[sujet+"_syn"][i]), df[sujet][j].extend(df[sujet+"_syn"][j]))
    return df_prox

def get_df_stats_prox(prox_mat_topic, prox_mat_operation):
    #prints mean min max std var for topics and operations of a set. you have to give the prox_dataframes as input
    #stats on the proximity score matrices
    val_top =  prox_mat_topic.values[np.triu_indices(len(prox_mat_topic), k = 1)] #extracting the values of the upper triangle of the matrix, diagonal not included
    val_op = prox_mat_operation.values[np.triu_indices(len(prox_mat_operation), k = 1)]

    #moyennes
    mean_top = np.mean(val_top)
    mean_op = np.mean(val_op)
    print(" moy topic = "+str(mean_top)+ ", moy operation = "+str(mean_op))

    #min-max
    min_top = np.min(val_top)
    min_op = np.min(val_op)
    print(" min topic = "+str(min_top)+ ", min operation = "+str(min_op))
    max_top = np.max(val_top)
    max_op = np.max(val_op)
    print(" max topic = "+str(max_top)+ ", max operation = "+str(max_op))


    #dispersion
    std_top = np.std(val_top)
    std_op = np.std(val_op)
    print(" std topic = "+str(std_top)+ ", std operation = "+str(std_op))

    var_top = np.var(val_top)
    var_op = np.var(val_op)
    print(" var topic = "+str(var_top)+ ", var operation = "+str(var_op))
    df = pd.DataFrame()
    df = df.append({"mean_top": mean_top,
        "min_top" : min_top,
        "max_top" : max_top,
        "std_top" : std_top,
        "var_top" : var_top,
        "mean_op": mean_op,
        "min_op" : min_op,
        "max_op" : max_op,
        "std_op" : std_op,
        "var_op" : var_op}, ignore_index=True)
    return df

























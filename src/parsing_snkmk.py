#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 08:51:09 2021

@author: marinedjaffardjy
"""
import numpy as np
import ntpath
import glob, os, shutil
import statistics
import snakemake.parser as snp
import get_characteristics as chr
import json
import requests
import re
from ratelimiter import RateLimiter
from difflib import SequenceMatcher
import search_biotools_dump as sb

###############################################################################
# Parsing the workflows with the snakemake parser

# using this path to test stuff out
# path = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/"
# file = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/davemcg/scEiaD/3.snakefile"
#print("hello")
# file2 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/bioxfu/RNA-Seq/1.snakefile"

#THIS FILE CONTAINS THE LIST OF THE EXCEPTIONS USED TO FILTER THE KEYWORDS
#TODO : import file and put it in a file
#TODO : enlever fonctions qui ne servent à rien

EXCEPTIONS = ["export", 'awk', 'sed', 'grep', "cmd", "module", "cat", "elif", "sort", "cd", "zcat",
              "rm", "for", "find", "java", "forgi", "sleep", "tabix", "zgrep", "wget", "mv", "mkdir", "echo",
              'FS', 'head', 'Rscript', "python", "jekyll","bgzip", "tr", "dot", "tRNA", "header", "fi", "then",
              "read", "do","else","cut", "wc", "tar", "gzip", "cool", "if", "turn", "git", "checkm", "cp", "make", "pour",
              "NR", "melt" , "read", "tail", "genes", "add", "bc", "scp", "scif", "uniq", "ln", "set", "zip", "time", "ls",
              "print","make", "pour", "source", "melt", "paste", "split", "layer", "touch" , "google-chrome", "query"
              "curl", "snps", "sh", "curl", "null", "join","config", "bin", "less", "comm", "vep", "which", "params",
              "next", "ml", "docker", "jupyter" , "date", "Date", "end", "END", "du", "Results", "idr", "latex", "tmp", "test",
              "command", "convert", "break", "sys", "root", "Root", "its" , "tmp", "nl", "ccs", "we", "We", "all", "ba",
              "mb", "lastal", "log", "must", "tl", "html", "all", "save", "tac", "use", "USE", "output", "PE", "output", "pre",
              "out", "file", "shell", "Track", "index", "In", "ref", "NF", "ION", "database", "ruby", "pos" , "load", "counts",
              "ra", "run", "are", "FROM", "nextflow", "DIP", "like", "to", "If", "chr", "env", "dt", "del", "clean", "fasta","sam",
              "my", "BF", "df", "view", "fq", "erp", "bp", "this", "The", "tag", "pip", "r1", "Dataset", "color" ]

# ajouté
def parse_simple(file):
    with ps.Snakefile(file, rulecount=0) as snakefile:
        automaton = ps.Python(snakefile)
    return automaton

def parse_wf(file):
    # parses a snakemake file using the snakemake parser
    # returns a dict of the number of rules and a string of the code for all the rules
    wf = snp.parse(file)
    wf_dict = {'rulecount': wf[2], 'compilation': wf[0], 'filename': file}
    return wf_dict

# à ajouter sous une autre forme

def parse_wf_folder(path):
    # explores a folder and its subdirectories to find all the snakefiles in it
    # returns a list of dicts of nb of rules and code for all the rules

    wf_files = []
    wf_names = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if (".snakefile" in file):  # this option can be changed depending on the snakefiles that you have, might change this later to put it in the params of the function
                wf_names.append(os.path.join(root, file))
    wf_names.sort()
    for file in wf_names:
        print(file)
        wf_files.append(parse_wf(os.path.join(root, file)))
    return wf_files

def parse_wf_folder_v2(path):
    # explores a folder and its subdirectories to find all the snakefiles in it
    # returns a list of dicts of nb of rules and code for all the rules

    wf_files = []
    wf_names = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if (".snakefile" in file):  # this option can be changed depending on the snakefiles that you have, might change this later to put it in the params of the function
                wf_names.append(os.path.join(root, file))
    wf_names.sort()
    for file in wf_names:
         #print(file)
        try :
            wf_files.append(parse_wf(os.path.join(root, file)))
        except SyntaxError :
            print("erreur de syntaxe sur le fichier "+str(os.path.join(root, file)))

    return wf_files
#ajout_erreurs
###############################################################################

# Getting additionnal information from the workflows

# aj
def get_rulelist(compilation):
    # separating the compilation data into a list of rules
    # returns the code of the list for each rule
    rules_list = compilation.split('@workflow.rule')[1:]
    return rules_list


# TODO : make another function for splitting the code of each rule into a dict for easier access of each part? not needed for now

#aj
def get_namelist(rules_list):
    # get the list of names from the rules list
    # return a list of names as a list of strings
    names = (rule[7:].split("'")[0] for rule in rules_list)
    return names

# a ajouter + tard
@RateLimiter(max_calls=5, period=1)
def get_info_biotools_tool(tool):
    # gets the informations on a given tool from the biotools API (for now considers the result with the best score as the most relevant == biotools' API chooses)
    # input : the name of the tool
    # returns a dict with the toolname, its biotools ID, collection ID, {topic}, {function}
    # TODO : make it so the toolname and the command have to be close enough otherwise nothing is returnes
    try:
        response = requests.get('https://bio.tools/api/tool/?format=json&name="' + tool + '"')
        if (response.status_code == 200):
            temp = json.loads(
                response.text)  # takes the first result returned by the biotools API, however can return stuff even if there is no real match with a tool
            # TODO : try to find the most relevant result
            # ISSUE WITH RScript and module for instance
            if (temp['count'] > 0):  # you can get a response 200 yet not have anything return by the biotools API
                temp = temp['list'][0]
                tool_info = {k: temp.get(k, None) for k in ('name', 'biotoolsID', 'topic', 'function')}
                # gets the toolnames, their biotools ID, their colletion ID, their {topic}, their {function} and gives the value none if the key doesn't exist (shoudldn't happen)
                # print('biotools ID ' + tool_info['biotoolsID'])
                return tool_info
            else:
                # print('empty list')
                return None
        else:
            # print ('status code 500')
            return None
    except TimeoutError:
        print('Max retry err')
        return None

# ne pas ajouter
def get_info_biotools_tool_v2(tool):
    # gets the informations on a given tool from the biotools API (for now considers the result with the best score as the most relevant == biotools' API chooses)
    # input : the name of the tool
    # returns a dict with the toolname, its biotools ID, colletion ID, {topic}, {function}
    response = requests.get('https://bio.tools/api/tool/?format=json&name=' + tool + '&sort=score&ord=desc')
    if (response.status_code == 200):
        temp = json.loads(response.text)  # takes the first result returned by the biotools API, however can return stuff even if there is no real match with a tool
        # TODO : try to find the most relevant result
        # ISSUE WITH RScript and module for instance
        if (temp['count'] > 0):  # you can get a response 200 yet not have anything return by the biotools API
            temp = temp['list'][0]
            tool_info = {k: temp.get(k, None) for k in ('name', 'biotoolsID', 'collectionID', 'topic', 'function')}
            # gets the toolnames, their biotools ID, their colletion ID, their {topic}, their {function} and gives the value none if the key doesn't exist (shoudldn't happen)
            # print('biotools ID ' + tool_info['biotoolsID'])
            return tool_info
        else:
            # print('empty list')
            return None
    else:
        # print ('status code 500')
        return None

# ajouté
def path_leaf(path):
    # extract last word from path
    # input : string path name
    # output string last element of path
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

# ajouté
def get_toolname_from_line(line):
    # get command from a line
    # input : a line (string)
    # output : a command word (string) or a tuple of two words (string, string) or none if the toolname is not relevant

    # TODO : solve issue in the case there are several keywords
    #print("line"+str(line))
    while(len(line)>1 and line[0]==' '):
        line = line[1:]
    line = line.replace('cmd ( "',"")
    tool = line.split(' ')[0].replace('\t', '').replace('"', '').replace("'", '')  # just for syntax, take out the tabs
    if ("=" in tool):
        tool.split('=')[0]
    if ("'" in tool):
        tool = tool.replace("'", '')
    # if the tool is in a command with a path ex /home/marine\.outil
    if ('/' in tool):
        tool = path_leaf(tool)
    if ('bowtie' in tool):
        tool = "bowtie"
    if ('htseq-count' in tool):
        tool = "htseqcount"
    #print(line)
    #print("tool "+tool)
     #look for second word :
    if len(line.split(' '))>1:
        sec_word = line.split(' ')[1] #weget the second word
        #if the second word starts with a letter and is bigger than one letter, we take it into account
        if(len(sec_word)>1):
            if (re.search('[a-zA-Z]', sec_word[0]) is not None and len(sec_word) > 1):
                #if the second word doesn't have an extension, we take it into account
                if ('.' not in sec_word) :
                    #print("second word " + sec_word)
                    return [tool,sec_word]
    return [tool]

# ajouté
def parse_lines(compil):
    # returns a list of lines from the shell data
    # input : shell content (string)
    # output : list of command lines as a list of strings
    lines_list = []
    for whole_line in compil.split('\n'):
        # check if the line has at least one character
        if (re.search('[a-zA-Z]', whole_line) is not None):  # checks that the line contains a character
            # split the line into several lines if there is a |, a ; or a >
            if any(ext in whole_line for ext in ['', '|', ';', '>']):
                whole_line = whole_line.split('|')
                for subline in whole_line:
                    subline = subline.split('>')
                    for subsubline in subline:
                        subsubline = subsubline.split(';')
                        for el in subsubline:
                            lines_list.append(el)
    #
    return lines_list

# ajouté
def get_toolnames(rulelist):
    # from the list of rules return a set of toolnames
    # TODO : make it a counter object so we have the info of the nb of times it was used
    # input : a list of string of the rules
    # returns a set of toolnames

    toolnames = []
    scripts_python = []
    scripts_R = []
    scripts_bash = []
    # print(rulelist)
    for rule in rulelist:
        #print(rule)
        if ('@workflow.shell' in rule):
            rule_lines = rule.split('@workflow.shell')[1]
            if ('@workflow.run' in rule_lines):
                rule_lines = rule_lines.split('@workflow.run')[0]
            # print("rule_lines "+rule_lines)
            lines_list = parse_lines(rule_lines)
            for line in lines_list:
                # print("line "+line)
                toolname = get_toolname_from_line(line)
                # print(toolname)
                if (re.search('[a-zA-Z]', toolname[0]) is not None and len(
                        toolname[0]) > 1):  # if the toolname is comprised of at least one letter and is of length sup to one

                    if (re.search('[a-zA-Z]', toolname[0][0]) is not None):  # if the toolname doesn't start with a special character
                        # # on cherche les scripts python
                        # if toolname[0] == "python":
                        #     scripts_python.append(toolname)
                        # # on cherche els scripts R
                        # if toolname [0]== "Rscript":
                        #     scripts_R.append(toolname)
                        # # on cherche les scripts bash
                        # if toolname[0][-3:]== ".sh":
                        #     scripts_bash.append(toolname)

                        if (toolname[0] not in EXCEPTIONS ):
                            if ('(' not in toolname[0] and '{' not in toolname[0] and '#' not in toolname[0] and "=" not in toolname[0] and '\\' not in toolname[0] and '.' not in toolname[0]):
                                toolnames.append(toolname)
                                # print("line " + line)

                                #print("toolname "+str(toolname))
    #        if ('shell ( "' in rule): #if we have a command line
    #            rule_lines = rule.split('shell ( "')[1]
    #            lines_list = parse_lines(rule_lines)
    #            for line in lines_list:
    #                print(line)
    #                toolname = get_toolname_from_line(line)
    #                print(toolname)
    #                toolnames.append(toolname)
    #return set(toolnames)
    #TODO : trouver un moyen de garantir l'unicité
    return toolnames
    # return toolnames, scripts_python, scripts_R, scripts_bash

# TODO : À AJOUTER PLUS TARD ?
def get_info_biotools_set_of_tools(set_tools):
    # queries the biotools api for a set of tools
    # intput : a set of strings of toolnames
    # output : a dict of dict of tools annotation, with the key being the given toolname
    dict_tools = {}
    for tool in set_tools:
        if (re.search('[a-zA-Z]',
                      tool) is not None):  # checks that the linetoolname contains a character, should be redundant but just in case
            if (tool not in ["module", "cat", "elif", "sort", "cd", "zcat", "rm", "zgrep", "wget", "mv", "mkdir",
                             "echo"]):
                # print('tool ' + tool)
                tool_info = get_info_biotools_tool(tool)
                if tool_info is not None:
                    # print('biotolsID ' + tool_info['biotoolsID'])
                    dict_tools.update({tool: tool_info})
    return dict_tools

# TODO : À AJOUTER DANS UN AUTRE ENDROIT - refaire avec le constr de tool ?
# TODO : TRIER EXCEPTIONS POUR LES AVOIR TOUTES AU MEME ENDROIT
def get_info_biotools_set_of_tools_dump(set_tools, treshold = 5):
    # queries the biotools dump for a set of tools
    #we can control the treshold of selectivity for the levensthein distance score
    # intput : a set of strings of toolnames
    # output : a dict of dict of tools annotation, with the key being the given toolname
    dict_tools = {}
    #TODO : faire une fct pour séparer le tri
    for tool in set_tools:
        if (re.search('[a-zA-Z]',tool[0]) is not None):  # checks that the linetoolname contains a character, should be redundant but just in case
            if (tool[0] not in ["module", "cat", "elif", "sort", "cd", "zcat", "rm", "zgrep", "wget", "mv", "mkdir",
                             "echo", "dot", "gunzip", "pandoc", "pdflatex", "python", "sleep", "done", "perl", "egrep", "tr", "rev", "jekyll", "rsync"]):
                # print('tool ' + str(tool))
                tool_info = sb.get_annotations_tool_dump(tool, treshold)
                if tool_info is not None:
                    # print('uri ' + tool_info['uri'])
                    if len(tool)==1:
                        dict_tools.update({tool[0]: tool_info})
                    else :
                        dict_tools.update({tool[0]+" "+tool[1]: tool_info})
    return dict_tools

# TODO : REFAIRE DS LE MAIN MAIS PR STOCKER UNE LISTE D'OBJETS TOOL ? ou d'identifiants ?
def attribute_info_biotools(wf, dict_tools):
    # from a dictionnary of the tools annotations add the annotations to the workflow
    # input a workflow, the dictionnary of toolnames
    # output the updated workflow with a list of tools annotations correspondong to its own tools
    list_toolinfo = []
    for tool in wf["extracted_toolnames"]:
        if len(tool) == 1:
            tool = tool[0]
        else:
            tool = tool[0]+" "+tool[1]
        if tool in dict_tools:
            list_toolinfo.append(dict_tools[tool])
    wf.update({"tools annotations": list_toolinfo})


def get_info_biotools(rules_list):
    # gets the informations on all the tools used in the workflow from the biotools API
    # input : the list of rules (list of strings)
    # returns a list of dict with the {toolnames, their biotools ID, their colletion ID, their {topic}, their {function}}

    # TODO : some additionnal metrics on the proportions of time that we have some success in querying biotools --
    tools = []
    for rule in rules_list:
        if ('@workflow.shell' in rule):  # if we have a command line
            # print('rule '+rule.split('shell ( "')[1])
            rule_lines = rule.split('@workflow.shell')[1]
            # TODO : faire un test case et une fonction pour "parser" les lignes
            for whole_line in rule_lines.split('\n'):
                # if there is a | in the line
                # TODO : essayer avec une ligne qui a ou ; ou >
                for line in whole_line.split('|'):
                    tool = get_toolname_from_line(line)
                    print(tool)
                    # if the toolname is not empty or "" or a parenthesis, query biotools
                    # making a manual exception for commands such as module and cat, rm, zgrep,wget, mv, mkdir
                    if (re.search('[a-zA-Z]', tool) is not None):  # checks that the linetoolname contains a character, should be redundant but just in case
                        if (tool not in ["module", "cat", "elif", "sort", "cd", "zcat", "rm", "zgrep", "wget", "mv",
                                         "mkdir", "echo"]):
                            #print('line ' + line)
                            #print('tool ' + tool)
                            tool_info = get_info_biotools_tool(tool)
                            if (tool_info != None):
                                tools.append(tool_info)
                                #print('appended ' + tool_info['biotoolsID'])
    # TODO : make a nonredundant list
    return tools

# TODO : à ajouter sous une autre forme ?
def get_info_biotools_dump(rules_list):
    # gets the informations on all the tools used in the workflow from the biotools API
    # input : the list of rules (list of strings)
    # returns a list of dict with the {toolnames, their biotools ID, their colletion ID, their {topic}, their {function}}

    # TODO : some additionnal metrics on the proportions of time that we have some success in querying biotools --
    tools = []
    for rule in rules_list:
        #print(rule)
        if ('@workflow.shell' in rule):  # if we have a command line
            #print('rule '+rule.split('shell ( "')[1])
            rule_lines = rule.split('@workflow.shell')[1]
            # TODO : faire un test case et une fonction pour "parser" les lignes
            for whole_line in rule_lines.split('\n'):
                # if there is a | in the line
                # TODO : essayer avec une ligne qui a ou ; ou >
                for line in whole_line.split('|'):
                    tool = get_toolname_from_line(line)
                    # if the toolname is not empty or "" or a parenthesis, query biotools
                    # making a manual exception for commands such as module and cat, rm, zgrep,wget, mv, mkdir
                    if (re.search('[a-zA-Z]',
                                  tool) is not None):  # checks that the linetoolname contains a character, should be redundant but just in case
                        if (tool not in ["module", "cat", "elif", "sort", "cd", "zcat", "rm", "zgrep", "wget", "mv",
                                         "mkdir", "echo"]):
                            #print('line ' + line)
                            #print('tool ' + tool)
                            tool_info = sb.get_annotations_tool_dump(tool)
                            if (tool_info != None):
                                tools.append(tool_info)
                                #print('appended ' + tool_info['biotoolsID'])
    # TODO : make a nonredundant list
    return tools


# Getting annotations information

# def bag_of_comments (rules_list):
# TODO : bag of comments

# Getting inputs/outputs names + extension
# kind of useless but could come in handy for annotation comparison

# ajouté
def get_inputnames(rules_list):
    # gets the names of the inputs from the rules list
    # returns a list of the names of the outputs as a string
    inputs = []
    for rule in rules_list:
        if ('@workflow.input(' in rule):
            inputs.append(rule.split('@workflow.input(')[1].split(')')[0])
    return inputs


# TODO : solve case of several inputs - same for outputs

def get_outputnames(rules_list):
    # gets the names of the inputs from the rules list
    # returns a list of the names of the outputs as a string
    outputs = []
    for rule in rules_list:
        if ('@workflow.output(' in rule):
            outputs.append(rule.split('@workflow.output(')[1].split(')')[0])
    return outputs


###############################################################################
# Creating the wf_informations 'object' and writing it into a file


def make_dict_wf(wf_inp):
    # makes a wf with the extracted toolnames
    wf_inp.update({'rules_list': get_rulelist(wf_inp['compilation'])})
    wf_inp.update({'rules_names': get_namelist(wf_inp['rules_list'])})
    wf_inp.update({'input_names': get_inputnames(wf_inp['rules_list'])})
    wf_inp.update({'output_names': get_outputnames(wf_inp['rules_list'])})
    wf_inp.update({'extracted_toolnames': get_toolnames(wf_inp['rules_list'])})
    # toolnames, scripts_python, scripts_R, bash = get_toolnames(wf_inp['rules_list'])
    # wf_inp.update({'extracted_toolnames': toolnames})
    # wf_inp.update({'scripts_python': scripts_python})
    # wf_inp.update({'scripts_R': scripts_R})
    # wf_inp.update({'scripts_bash': bash})
    return wf_inp

#USING A LOCAL DUMP

def make_dict_wf_biotools_dump_wf_by_wf(wf_inp):
    # from the initial dict made using the snakemake parser
    # returns a dict with : - filename : a string of the name of the snakefile with its path
    #                      - rulecount : tan int of he number of rules
    #                      - compilation : a string of the code for all the rules
    #                      - rules_list : a list of strings for the code of all the rules (split of the compilation above)
    #                      - rule_names : a list of strings for all the names of the rules
    #                      - input_names : a list of strings for all the output names
    #                      - output_names : a list of strings for all the input names
    #                      - tools : list of dict with the {toolnames, their biotools ID, their colletion ID, their {topic}, their {function}}

    wf_inp.update({'rules_list': get_rulelist(wf_inp['compilation'])})
    wf_inp.update({'rules_names': get_namelist(wf_inp['rules_list'])})
    wf_inp.update({'input_names': get_inputnames(wf_inp['rules_list'])})
    wf_inp.update({'output_names': get_outputnames(wf_inp['rules_list'])})
    wf_inp.update({'tools_list': get_info_biotools_dump(wf_inp['rules_list'])})
    return wf_inp


def make_dict_wf_with_biotools_dump(wf_tab):
    # for a list of wf, updates the list of wf to add the keys rules list, rule names, unput names, output names
    # also takes the extracted toolnames set
    total_extracted_tools = []
    for wf in wf_tab:
        wf = make_dict_wf(wf)
        total_extracted_tools.extend(wf['extracted_toolnames'])
    total_extracted_tools = set(tuple(i) for i in total_extracted_tools)
    dict_tools = get_info_biotools_set_of_tools_dump(total_extracted_tools)
    #print(dict_tools)
    for wf in wf_tab:
        attribute_info_biotools(wf, dict_tools)
    return dict_tools, total_extracted_tools

def make_dict_wf_with_biotools_dump(wf_tab):
    # for a list of wf, updates the list of wf to add the keys rules list, rule names, unput names, output names
    # also takes the extracted toolnames set
    total_extracted_tools = []
    for wf in wf_tab:
        wf = make_dict_wf(wf)
        total_extracted_tools.extend(wf['extracted_toolnames'])
    total_extracted_tools = set(tuple(i) for i in total_extracted_tools)
    dict_tools = get_info_biotools_set_of_tools_dump(total_extracted_tools)
    #print(dict_tools)
    for wf in wf_tab:
        attribute_info_biotools(wf, dict_tools)
    return dict_tools, total_extracted_tools

# USING THE BIOTOOLS API
def make_dict_wf_biotools_by_wf(wf_inp):
    # from the initial dict made using the snakemake parser
    # returns a dict with : - filename : a string of the name of the snakefile with its path
    #                      - rulecount : tan int of he number of rules
    #                      - compilation : a string of the code for all the rules
    #                      - rules_list : a list of strings for the code of all the rules (split of the compilation above)
    #                      - rule_names : a list of strings for all the names of the rules
    #                      - input_names : a list of strings for all the output names
    #                      - output_names : a list of strings for all the input names
    #                      - tools : list of dict with the {toolnames, their biotools ID, their colletion ID, their {topic}, their {function}}

    wf_inp.update({'rules_list': get_rulelist(wf_inp['compilation'])})
    wf_inp.update({'rules_names': get_namelist(wf_inp['rules_list'])})
    wf_inp.update({'input_names': get_inputnames(wf_inp['rules_list'])})
    wf_inp.update({'output_names': get_outputnames(wf_inp['rules_list'])})
    wf_inp.update({'tools_list': get_info_biotools(wf_inp['rules_list'])})
    return wf_inp


def make_dict_wf_with_biotools(wf_tab):
    # for a list of wf, updates the list of wf to add the keys rules list, rule names, unput names, output names
    # also takes the extracted toolnames set
    total_extracted_tools = set()
    for wf in wf_tab:
        wf = make_dict_wf(wf)
        total_extracted_tools.update(wf['extracted_toolnames'])
    #get_unique tool list
    dict_tools = get_info_biotools_set_of_tools(total_extracted_tools)
    for wf in wf_tab:
        attribute_info_biotools(wf, dict_tools)
    return dict_tools, total_extracted_tools


# Write the json files in txt files (in order to not redo everything everytime)
def write_json_in_txt(wf_dict, outputdir, inputdir):
    # from a wf dict with all the info write everything in a text file
    with open(outputdir + wf_dict['filename'].split(inputdir)[1] + '.txt', 'w') as outfile:
        json.dump(outfile)




# if __name__ == '__main__':
#     path_800 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/data/"
#
#     # workflow parsing
#     wf_800 = parse_wf_folder(path_800)
#
#     dict_tool, total_extr = ps.make_dict_wf_with_biotools_dump(wf_800)
#
#     write_json_in_txt(wf_dict, outputdir, inputdir)
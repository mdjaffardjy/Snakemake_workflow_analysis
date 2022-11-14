#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 08:51:09 2021

@author: marinedjaffardjy
"""

import ntpath
import os
import re
import Snakemake.parser_modif as snp
import search_biotools_dump as sb


# The exceptions list for finding tools

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



###########################################
## FINDING AND PARSING THE SNAKEFILES
###########################################

def parse_wf(file):
    # parses a snakemake file using the snakemake parser
    # returns a dict of the number of rules and a string of the code for all the rules
    wf = snp.parse(file)
    wf_dict = {'rulecount': wf[2], 'compilation': wf[0], 'filename': file}
    return wf_dict


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
         #print(file)
        try :
            wf_files.append(parse_wf(os.path.join(root, file)))
        except :
            print("parsing error on file "+str(os.path.join(root, file)))

    return wf_files

###########################################
## ANALYZING THE BASIC PARSING INFORMATION
###########################################

def get_rulelist(compilation):
    # separating the compilation data into a list of rules
    # returns the code of the list for each rule
    rules_list = compilation.split('@workflow.rule')[1:]
    return rules_list



def get_name(rule):
    # get the names of a rule
    # return a names as a strings
    name = rule[7:].split("'")[0]
    return name


def get_inputnames(rule):
    # get the inputs of a rule
    # return inputs as a string
    
    inputs = ""
    if ('@workflow.input(' in rule):
        inputs = rule.split('@workflow.input(')[1].split(')')[0]
    return inputs



def get_outputnames(rule):
    # get the outputs of a rule
    # return outputs as a string
    
    outputs = ""
    if ('@workflow.output(' in rule):
        outputs = rule.split('@workflow.output(')[1].split(')')[0]
    return outputs


###########################################
## FINDING THE TOOLS IN THE RULES
###########################################

def path_leaf(path):
    # extract last word from path
    # input : string path name
    # output string last element of path
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def get_toolname_from_line(line):
    # get command from a line
    # input : a line (string)
    # output : a command word (string) or a tuple of two words (string, string) or none if the toolname is not relevant

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
    # look for second word :
    if len(line.split(' '))>1:
        sec_word = line.split(' ')[1] #weget the second word
        #if the second word starts with a letter and is longer than one letter, we take it into account
        if(len(sec_word)>1):
            if (re.search('[a-zA-Z]', sec_word[0]) is not None and len(sec_word) > 1):
                #if the second word doesn't have an extension, we take it into account
                if ('.' not in sec_word) :
                    #print("second word " + sec_word)
                    return [tool,sec_word]
    return [tool]


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
    
    return lines_list


def get_toolnames(rulelist):
    # from the list of rules return a set of toolnames
    # input : a list of string of the rules
    # returns a set of toolnames

    toolnames = []
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
                        if (toolname[0] not in EXCEPTIONS ):
                            if ('(' not in toolname[0] and '{' not in toolname[0] and '#' not in toolname[0] and "=" not in toolname[0] and '\\' not in toolname[0] and '.' not in toolname[0]):
                                toolnames.append(toolname)
                                
    return toolnames

###########################################
## SEARCHING FOR THE TOOLS IN THE BIOTOOLS DUMP
###########################################

def get_info_biotools_set_of_tools_dump(set_tools, treshold = 5):
    # queries the biotools dump for a set of tools
    #we can control the treshold of selectivity for the levensthein distance score
    # intput : a set of strings of toolnames
    # output : a dict of dict of tools annotation, with the key being the given toolname
    dict_tools = {}
    for tool in set_tools:
        if (re.search('[a-zA-Z]',tool[0]) is not None):  # checks that the linetoolname contains a character, should be redundant but just in case
            if (tool[0] not in ["module", "cat", "elif", "sort", "cd", "zcat", "rm", "zgrep", "wget", "mv", "mkdir",
                             "echo", "dot", "gunzip", "pandoc", "pdflatex", "python", "sleep", "done", "perl", "egrep", "tr", "rev", "jekyll", "rsync"]):
                tool_info = sb.get_annotations_tool_dump(tool, treshold)
                if tool_info is not None:
                    if len(tool)==1:
                        dict_tools.update({tool[0]: tool_info})
                    else :
                        dict_tools.update({tool[0]+" "+tool[1]: tool_info})
    return dict_tools

def get_info_biotools_dump(rules_list):
    # gets the informations on all the tools used in a rule from the biotools API
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
    return tools

###########################################
## CREATING WORFLOW DICTIONNARIES WITH EXTRACTED INFORMATION
###########################################
  

def make_dict_wf_biotools_dump_wf_by_wf(wf_inp):
    # from the initial dict made using the snakemake parser for a single workflow
    # search biotools workflow per workflow
    #                      - filename : a string of the name of the snakefile with its path
    #                      - rulecount : tan int of he number of rules
    #                      - compilation : a string of the code for all the rules
    #                      - rules_list : a list of dict for the rules with :
    #                               - rule name, 
    #                               - rule codes,  
    #                               - inputs,  
    #                               - outputs   
    #                               - tools : list of dict with : 
    #                                   - toolnames, their biotools ID,  
    #                                   - biotools ID,  
    #                                   - topics,  
    #                                   - their function (operations, inputs and outputs)
    
    rules_codes = get_rulelist(wf_inp['compilation'])
    rules_list_dict = []
    for rule in rules_codes :
        name_rule = get_name(rule)
        input_rule = get_inputnames(rule)        
        output_rule = get_outputnames(rule)
        rule_tools = get_info_biotools_dump([rule])
        rules_list_dict.append({"rule_name" : name_rule,
                                "rule_code" : rule,
                                "input_rule" : input_rule,
                                "output_rule" : output_rule,
                                "tools": rule_tools})
    wf_inp.update({"rules" : rules_list_dict})
    return wf_inp



def make_dict_wf(wf_inp):
    # makes a wf dictionnary with all 'basic attributed" as well extracted toolnames (candidates toolnames to be matched against biotools)
    rules_codes = get_rulelist(wf_inp['compilation'])
    rules_list_dict = []
    for rule in rules_codes :
        name_rule = get_name(rule)
        input_rule = get_inputnames(rule)        
        output_rule = get_outputnames(rule)
        extracted_toolnames = get_toolnames([rule])
        rules_list_dict.append({"rule_name" : name_rule,
                                "rule_code" : rule,
                                "input_rule" : input_rule,
                                "output_rule" : output_rule,
                                "extracted_toolnames": extracted_toolnames})
    wf_inp.update({"rules" : rules_list_dict})
    return wf_inp

def attribute_info_biotools(wf, dict_tools):
    # from a dictionnary of the tools annotations add the annotations to the workflow
    # input a workflow, the dictionnary of toolnames
    # output the updated workflow with a list of tools annotations correspondong to its own tools
    list_toolinfo = []
    for rule in wf["rules"] :
        for tool in rule["extracted_toolnames"]:
            if len(tool) == 1:
                tool = tool[0]
            else:
                tool = tool[0]+" "+tool[1]
            if tool in dict_tools:
                list_toolinfo.append(dict_tools[tool])
                rule.update({"tools annotations": list_toolinfo})
    



def make_dict_wf_with_biotools_dump(wf_tab):
    # for a list of wf, updates the list of wf to add the keys rules list, rule names, unput names, output names
    # same as make_dict_wf_biotools_dump_wf_by_wf but is better to use when dealing with a large amount of workflow as it queries the dump only once for all wf and not once for every workflow
    
    # modifies the workflow dictionnary to a dictionnary with the keys :
    #                      - filename : a string of the name of the snakefile with its path
    #                      - rulecount : tan int of he number of rules
    #                      - compilation : a string of the code for all the rules
    #                      - rules_list : a list of dict for the rules with :
    #                               - rule name, 
    #                               - rule codes,  
    #                               - inputs,  
    #                               - outputs   
    #                               - tools : list of dict with : 
    #                                   - toolnames, their biotools ID,  
    #                                   - biotools ID,  
    #                                   - topics,  
    #                                   - their function (operations, inputs and outputs)

    total_extracted_tools = []
    
    
    for wf in wf_tab:
        print(wf["filename"])
        wf = make_dict_wf(wf)
        for rule in wf["rules"]:
            total_extracted_tools.extend(rule['extracted_toolnames'])
    # match the extracted toolnames of all workflows to biotools dump
    total_extracted_tools = set(tuple(i) for i in total_extracted_tools)
    dict_tools = get_info_biotools_set_of_tools_dump(total_extracted_tools)
    #print(dict_tools)
    #attribute the tool info to all workflows
    for wf in wf_tab:
        attribute_info_biotools(wf, dict_tools)
    return dict_tools, total_extracted_tools







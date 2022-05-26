#!/usr/bin/env python

# Copyright (C) 2019 Thermo Fisher Scientific. All Rights Reserved

import glob
import sys
import subprocess
import json
import os
import re
import shutil
from ion.plugin import *
from django.utils.functional import cached_property
from django.conf import settings
from django.template.loader import render_to_string

'''
change log
  2022-5-22: 1.0.0.0

  2022-5-25: 1.0.0.1
    1)modify header;
    2)add some cols:
      * P1/P2_Ratio
      * Loading(%)
      * Enrichment(%)
      * Polyclonal(%)
      * Low Quality(%)
      * Primer Dimer(%)
  
  2022-5-26: 1.0.0.2
    1) add genomeCoverage(%) col [1-consN]
  ...
'''

class ReportSummary_2019nCoV(IonPlugin):
  version = '1.0.0.2'
  author  = "longfei.fu@thermofisher.com"
  date    = "2022-5-12"
  runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]

  # a simple cached version of the start plugin property
  @cached_property
  def startplugin_json(self):
    return self.startplugin

  def launch(self,data=None):
    net_location = self.startplugin_json['runinfo']['net_location'] # http://9RPWLN2
    plugin_result_dir = self.startplugin_json['runinfo']['plugin'].get('results_dir') # /results/analysis/output/Home/S5yanzheng-20211223-chip2-MeanAccuracy_v2_482/plugin_out/variantCaller58WES_out.1783
    report_dir = self.startplugin_json['runinfo']['report_root_dir'] # /results/analysis/output/Home/S5yanzheng-20211223-chip2-MeanAccuracy_v2_482
    abs_path = os.path.abspath(__file__)
    this_dir = os.path.dirname(abs_path)
    
    cmd = "python %s/ReportSummary_2019nCoV_plugin.py -report_dir %s -outdir %s" % (this_dir,report_dir,plugin_result_dir)
    print "cmd is: %s" % (cmd)
    os.system(cmd)

    url_root = self.startplugin_json['runinfo']['url_root'] # /output/Home/S5yanzheng-20211223-chip2-MeanAccuracy_v2_482
    file_path = "%s/plugin_out/%s" % (url_root,os.path.basename(plugin_result_dir))
    print(file_path)

    with open("ReportSummary_2019nCoV_block.html","w") as f:
      # first print which plugin was used
      plugin_log = "%s/plugin.log" % (plugin_result_dir)
      print(plugin_log)
      of_p = open(plugin_log,'r')
      lines = of_p.readlines()
      f.write("<html><body>\n")
      f.write("<ul>\n")
      for line in lines:
        v = '  <li style="color:Olive">%s</li>' % (line.rstrip())
        f.write(v+'\n')
      f.write("</ul>\n")
         
      for file in glob.glob('*.summary.xls'):
        print(file)
        f.write('<a href="%s">%s</a><br>\n'
              % (os.path.join(net_location,file_path,file),file))
        f.write('</body></html>')
    
    return True


if __name__ == "__main__":
    PluginCLI()

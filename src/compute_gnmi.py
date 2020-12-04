#!/usr/bin/env python

# Copyright 2020 MaaT Pharma
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
import subprocess  
import pandas as pd
import re
import numpy as np

def compute_gnmi_matrix(cnl_files,cnl_files_label, tool_path):

	NMImax_results = np.empty((len(cnl_files_label),len(cnl_files_label)))
	NMImax_results[:] = np.nan

	for cnlF in range(len(cnl_files)):
		for cnlF_comp in range(len(cnl_files)) :

			cmd = "{} {} {} | grep -v 'WARNING'".format(os.path.join(tool_path,"gecmi"),cnl_files[cnlF], cnl_files[cnlF_comp])

			# print(cmd)
			try: 
				output = subprocess.check_output(cmd,shell=True)
 
				NMI_find_max = re.match('^(.*)$',(output.decode('ascii')).split("\n")[0])
				if NMI_find_max and len(NMI_find_max.group()) : 
					NMImax_results[cnlF,cnlF_comp] = NMI_find_max.group(1)
				else : 
					NMImax_results[cnlF,cnlF_comp] = float('nan')

			except subprocess.CalledProcessError as e:
				print("WARNING: error computing GNMI between {} and {}".format(cnl_files_label[cnlF], cnl_files_label[cnlF_comp])) 
				
				NMImax_results[cnlF,cnlF_comp] = float('nan')

	df_NMImax_results= pd.DataFrame(NMImax_results,columns=cnl_files_label, index=cnl_files_label)

	return df_NMImax_results


def print_gnmi_results(df_NMImax_results, output_dir) : 
	df_NMImax_results.to_csv(os.path.join(output_dir, "GNMI_max.tsv"), sep='\t')
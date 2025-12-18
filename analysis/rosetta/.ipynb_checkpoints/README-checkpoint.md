ddG predictions for GFP proteins were calculated using **Rosetta** software, with the following command (replace "/rosetta/rosetta-3.10/" with your Rosetta location):

`/rosetta/rosetta-3.10/main/source/bin/ddg_monomer.linuxgccrelease  @flags.txt -in:file:extra_res_fa GYS.params -ddg::mut_file mutfile_0.txt -ddg::out all_out.txt`
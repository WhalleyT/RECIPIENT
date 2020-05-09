# Bacterial Vaccine Design

### Download
The tool can be downloaded by running `git clone https://github.com/whalleyt/bacterial_vaccine_design`. This directory contains the bare bones to run the tool within docker containers.

The Docker containers can be downloaded as follows:
`docker pull twhalley93/loctree3:latest`
`docker pull twhalley93/bvd:latest`

The Dockerfiles, along with the source code of the entire script (should you wish to run the code locally) can be found [here](github.com/whalleyt/bvd_source). **please** consider using the Dockerfiles provided, this is because a number of tools have to be manually downloaded manually due to licensing constraints.

### Run
The tool can be ran using nextflow as follows:
`nextflow pipeline.nf --fasta_folder <folder> --outdir <outdir>`

Help can be accessed by `nextflow pipeline -h`.

### License
All tools intgegrated in this software, and the source code are published under the GNU GPL; with the exception of NetMHCIpan and NetMHCIIpan. Please find the academic license in the NETMHC_LICENSE file. If you are a non-academic user, please contact the authors of NetMHC before proceeding. They can be found at: software@cbs.dtu.dk.

The license documentation can be found in LICENSE.

### Contact
please feel free to contact us at Whalleyt@cardiff.ac.uk, alternatively make a pull request.

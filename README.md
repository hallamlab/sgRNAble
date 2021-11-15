# sgRNAble

__sgRNAble__ is a tool for high-throughput design of sgRNA libraries targeting selected genes or whole genomes, while considering both on-target binding potential and off-target effects of a given sgRNA in a user-defined genome.

## Basic Usage

```
 sgrnable -t TARGET_SEQUENCE.fasta -g GENOME_SEQUENCE.fasta
          [GENOME_SEQUENCE.fasta ...] [-a AZIMUTH_CUTOFF]
          [-c COPY_NUMBER [COPY_NUMBER ...]]
          [-o OUTPUT_DIR] [-p PURPOSE]
          [-th NUM_THREADS] [-m MAX_MEMORY] [-v]
```

## Documentation

All documentation for the current release is accessible through the [sgRNAble wiki](https://github.com/hallamlab/sgRNAble/wiki).

## Authors
* [Siddarth Raghuvanshi](https://github.com/Siddarth-Raghuvanshi)
* [Ahmed Abdelmoneim](https://github.com/AhmedAbdelmoneim)
* [Avery Noonan](https://github.com/Noonanav)

## Contact

Need something? Please contact Professor Steven Hallam and Avery Noonan at: shallam@mail.ubc.ca and avery.noonan@ubc.ca

## References

Farasat, I., & Salis, H. M. (2016). A Biophysical Model of CRISPR/Cas9 Activity for Rational Design of Genome Editing and Gene Regulation. _PLOS Computational Biology_, __12__(1), e1004724. doi:10.1371/journal.pcbi.1004724

Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan, K. F., . . . Root, D. E. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. _Nature biotechnology_, __34__(184). doi:10.1038/nbt.3437

qiime_pipe
==========

Ruby script to execute all commands in the QIIME overview tutorial, but including denoising and chimera check.

Works with QIIME 1.7

```
Usage: qiime_pipe.rb [options]
    -h, --help                       Display this screen
    -s, --file_sff <file>            SFF file to process
    -m, --file_map <file>            Mapping file to process
    -o, --dir_out <dir>              Output directory
    -f, --force                      Force overwrite log file and output directory
    -d, --denoise                    Denoise data
    -c, --chimera                    Chimera filter data
    -D, --chimera_db <file>          Chimere database (/home/maasha/Install/QIIME1.6/data/Gold/gold.fa)
    -b, --barcode_size <int>         Size of barcodes used (10)
    -C, --cpus <int>                 Number of CPUs to use (1)
    -e, --email <string>             Send email alert
```

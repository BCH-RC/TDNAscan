# Installing TDNAscan 

## Operating System Requirements

TDNAscan has been tested on the following Linux distributions:

* Ubuntu 14.04 LTS
* Ubuntu 16.04 LTS
* CentOS 7.3
* Debian 7 "Wheezy"
* Debian 8 "Jessie"

## Install from Source

1. Run: `git clone` (or download tdnascan.tar.gz from Releases and run `tar -xzvf fnbtools.tar.gz`)
2. Run: `sudo ./Install.sh`
3. If prompted with "Would you like to configure as much as possible automatically? [yes]", type **yes** and then press Enter. This will automatically configure Perl's CPAN utility so that additional Perl modules can be installed.
4. If prompted with "Would you like me to automatically choose some CPAN mirror sites for you? [yes]", type **yes** and then press Enter.  This will automatically configure Perl's CPAN utility with a mirror site from which it can download Perl modules.

**Optional** - If you would like to enable visualization of the results, continue with steps 5 - 7:

5. Run: `echo 'export PATH=/usr/local/circos/current/bin:$PATH' >> ~/.bashrc`
6. Run: `. ~/.bashrc`
7. Run `circos -modules` to verify that all Perl modules that are required to run Circos for visualizations are installed. If any are reported as missing, you must install them before attempting to visualize FNBTools results.
    * A helper script called `ReinstallCircosPerlModules.sh` has been provided to assist with installing Perl modules that fail to install during the main installation.
    * Run: `sudo ./ReinstallCircosPerlModules.sh` to re-attempt to install missing Perl modules that Circos requires.



# Using TDNAscan 

# Contact

* Dr. Liang Sun    lsun@noble.org
* Yinbing    ybge@noble.org

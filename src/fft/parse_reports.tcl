# parse_reports.tcl --
#
###Author
# Hugues ALMORIN (hugues.almorin@arelis.com)
###Abstract
# This file, allows to create fft run reports for all selected configurations
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3
#

source "../../scripts/report_parser.tcl"

proc build_csynth_table { } {
    set hdl_lang "vhdl"
    set fileout "./report.txt"
    set path "./hls_runs"
    set arch "fft_sbg_radix"
    set top "fft"
    set targets_name [list "Kintex_U+"]
    #set quantiz {"8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24"}
    set quantiz {"18"}
    #set tw_width {"10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26"}
    set tw_width {"20"}
    set freq {"50MHz" "100MHz" "200MHz" "300MHz"}
    set int_freq {50 100 200 300}
    set radix {2 4 8}
    # set radix {4 8}
    set NFFT { {"64" "128" "256" "512" "1024" "2048" "4096"} {"64" "256" "1024" "4096"} {"64" "512" "4096"} }
    # set NFFT { {"64" "256" "1024" "4096"} {"64" "512" "4096"} }
    set SW { {2 4 8} {4 8} {8} }
    # set SW { {4 8} {8} }
    set GROUP { {"1" "2" "3" "4" "5" "6"} {"1" "2" "3"} {"1" "2"} }
    # set GROUP { {"1" "2" "3"} {"1" "2"} }
    set ftime [clock format [clock seconds] -format "%H.%M.%S-%d.%m.%Y"]
    set csynth_fileout "./${arch}_csynth_report_${ftime}.txt"
    set export_fileout "./${arch}_export_report_${ftime}.txt"
    set path "./hls_runs"

    set arch_param_names [list "Arch" "Target" "NFFT" "Freq" "Quantif" "SW" "Group" "Radix"]
    set arch_param_names_format "%20s %12s %6s %6s %8s %6s %6s %6s"

    set csynth_fd_out [open ${csynth_fileout} w]
    set export_fd_out [open ${export_fileout} w]

    puts "------------------------------"
    puts "CSYNTH REPORT"
    puts "------------------------------"
    report_parser::build_csynth_header ${csynth_fd_out} ${arch_param_names} ${arch_param_names_format}
    foreach r $radix nfft $NFFT sw $SW group $GROUP {
        foreach g $group {
            foreach s $sw {
                foreach t $targets_name {
                    foreach n $nfft {
                        foreach f $freq intf $int_freq {
                            foreach q $quantiz {
                                report_parser::build_csynth_line ${csynth_fd_out} "${path}/${arch}_N${n}_Q2_${q}_SW${s}_G${g}_R${r}/${t}_${f}/syn/report/csynth.xml" [list $arch $t $n $intf ${q} ${s} ${g} ${r}] $arch_param_names ${arch_param_names_format}
                            }
                        }
                    }
                }
            }
        }
    }

    puts "\n------------------------------"
    puts "EXPORT REPORT"
    puts "------------------------------"
    report_parser::build_export_header ${export_fd_out} ${arch_param_names} ${arch_param_names_format}
    foreach r $radix nfft $NFFT sw $SW group $GROUP {
        foreach g $group {
            foreach s $sw {
                foreach t $targets_name {
                    foreach n $nfft {
                        foreach f $freq intf $int_freq {
                            foreach q $quantiz {
                                report_parser::build_export_line ${export_fd_out} "${path}/${arch}_N${n}_Q2_${q}_SW${s}_G${g}_R${r}/${t}_${f}/impl/report/${hdl_lang}/${top}_export.xml" [list $arch $t $n $intf ${q} ${s} ${g} ${r}] $arch_param_names ${arch_param_names_format}
                            }
                        }
                    }
                }
            }
        }
    }
    
    close $csynth_fd_out
    close $export_fd_out
}

build_csynth_table
quit
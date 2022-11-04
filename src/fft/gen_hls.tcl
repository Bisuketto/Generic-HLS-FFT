# gen_hls.tcl --
#
###Author
# Hugues ALMORIN (hugues.almorin@arelis.com)
###Abstract
# This file, called within Xilinx Vitis HLS 2020.2, allows to run the fft for all selected configurations
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3
#


source "../../scripts/report_parser.tcl"

set path "./hls_runs"
file mkdir $path
cd $path

set hdl_lang "vhdl"
set arch "fft_sbg_radix"
set top "fft"
set targets_name [list "Kintex_U+"]
set targets_id [list xcku5p-sfvb784-2-i]
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

foreach r $radix nfft $NFFT sw $SW group $GROUP {
    foreach g $group {
        foreach s $sw {
            foreach n $nfft {
                foreach q $quantiz tw $tw_width {
                    open_project -reset ${arch}_N${n}_Q2_${q}_SW${s}_G${g}_R${r}
                    set_top ${top}

                    set flags "-DNFFT=${n} \
                                -DCOS_LUT=COS_LUT_N${n}_Q2_${q} \
                                -DSIN_LUT=SIN_LUT_N${n}\_Q2_${q} \
                                -DDIGREV_CONFIG=DIG_REV_PERM_CONFIG_N${n}_SW${s}_R${r} \
                                -DDIGREV_NUMSTAGE=DIGIT_REV_NUM_STAGE_N${n}_SW${s}_R${r} \
                                -DSTRIDE_CONFIG=STRIDE_PERM_CONFIG_N${n}_SW${s}_R${r} \
                                -DSTRIDE_NUMSTAGE=STRIDE_PERM_SWITCH_NUM_STAGE_N${n}_SW${s}_R${r} \
                                -DSTREAMING_WIDTH=${s} \
                                -DTW_WIDTH=${tw} \
                                -DGROUP_SIZE=${g} \
                                -DCORE_RADIX=${r}"

                    add_files ../module/fft_sbg_radix_hlstop.cpp -cflags "${flags}" -csimflags "${flags}"
                    
                    add_files -tb ../tb/fft_sbg_radix_tb.cpp -cflags "${flags}" -csimflags "${flags}"
                    add_files -tb ../../common/fft_fftw3.cpp -cflags "-lfftw3 ${flags}" -csimflags "-lfftw3 ${flags}"

                    foreach f $freq {
                        foreach tname $targets_name tid $targets_id {
                            open_solution -reset "${tname}_${f}" -flow_target vivado
                            set_part $tid
                            create_clock -period ${f}

                            csim_design
                            csynth_design
                            cosim_design -rtl ${hdl_lang}
                            export_design -format ip_catalog -flow impl -rtl ${hdl_lang}

                            close_solution
                        }
                    }

                    close_project
                }
            }
        }
    }
}

proc build_csynth_table { } {
    variable hdl_lang
    variable arch
    variable targets_name
    variable NFFT
    variable freq
    variable int_freq
    variable quantiz
    variable SW
    variable GROUP
    variable radix
    variable top
    variable csynth_fileout
    variable export_fileout

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
                                report_parser::build_csynth_line ${csynth_fd_out} "./${arch}_N${n}_Q2_${q}_SW${s}_G${g}_R${r}/${t}_${f}/syn/report/csynth.xml" [list $arch $t $n $intf ${q} ${s} ${g} ${r}] $arch_param_names ${arch_param_names_format}
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
                                report_parser::build_export_line ${export_fd_out} "./${arch}_N${n}_Q2_${q}_SW${s}_G${g}_R${r}/${t}_${f}/impl/report/${hdl_lang}/${top}_export.xml" [list $arch $t $n $intf ${q} ${s} ${g} ${r}] $arch_param_names ${arch_param_names_format}
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

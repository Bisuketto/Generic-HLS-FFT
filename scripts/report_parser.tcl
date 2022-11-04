# gen_hls.tcl --
#
###Author
# Hugues ALMORIN (hugues.almorin@arelis.com)
###Abstract
# This file, allows to extract all run reports from a Xilinx Vitis HLS 2020.2 solution
#
# This file is released under the GNU GENERAL PUBLIC LICENSE Version 3
#

package require tdom

namespace eval report_parser {
    proc get_fields { fields node } {
        foreach field $fields {
            set fnode [$node selectNodes $field]
            if {[llength $fnode] == 0} {
                lappend values "N/A"
            } else {
                set fvalue [$fnode text]
                if {$fvalue == ""} {
                    lappend values "N/A"
                } else {
                    lappend values $fvalue
                }
            }
        }
        return $values
    }

    set resources_path "AreaEstimates/Resources"
    set required_resources [list "DSP" "BRAM_18K" "LUT" "FF" "URAM"]
    set required_resources_format "%8s %8s %8s %8s %8s"
    
    set timing_path "PerformanceEstimates/SummaryOfTimingAnalysis"
    set required_timing [list "EstimatedClockPeriod"]
    set required_timing_format "%20s"

    set latency_path "PerformanceEstimates/SummaryOfOverallLatency"
    set required_latency [list "Best-caseLatency" "Average-caseLatency" "Worst-caseLatency" "Interval-min" "Interval-max"]
    set required_latency_format "%16s %19s %17s %12s %12s"

    proc build_csynth_header { fd_out arch_param_names arch_param_names_format } {

        variable required_resources
        variable required_resources_format
        
        variable required_timing
        variable required_timing_format

        variable required_latency
        variable required_latency_format

        set line_format "${arch_param_names_format} ${required_resources_format} ${required_timing_format} ${required_latency_format}"
        set table_header [concat $arch_param_names $required_resources $required_timing $required_latency]
        
        puts [format $line_format {*}$table_header]
        puts $fd_out [format $line_format {*}$table_header]
    }

    proc build_csynth_line { fd_out fname arch_param arch_param_names arch_param_names_format} {
        variable resources_path
        variable required_resources
        variable required_resources_format
        
        variable timing_path
        variable required_timing
        variable required_timing_format

        variable latency_path
        variable required_latency
        variable required_latency_format

        set line_format "${arch_param_names_format} ${required_resources_format} ${required_timing_format} ${required_latency_format}"

        if { [catch {open $fname} xmlfd] } {
            set nparams [expr [llength $required_resources] + [llength $required_timing] + [llength $required_latency]]
            for {set i 1} {$i <= $nparams} {incr i} {
                append results " FAILED"
            }
            set line [concat $arch_param $results]
            puts [format $line_format {*}$line]
            puts $fd_out [format $line_format {*}$line]
        } else {        
            set xmlcontent [read $xmlfd]
            close $xmlfd

            set xmldoc [dom parse $xmlcontent]
            set root [$xmldoc documentElement]

            set resources_node [$root selectNodes $resources_path]
            set timing_node [$root selectNodes $timing_path]
            set latency_node [$root selectNodes $latency_path]

            set line [concat $arch_param                            \
                [get_fields $required_resources $resources_node]    \
                [get_fields $required_timing $timing_node]          \
                [get_fields $required_latency $latency_node]        ]
            puts [format $line_format {*}$line]
            puts $fd_out [format $line_format {*}$line]
        }
    }

    set impl_resources_path "AreaReport/Resources"
    set required_impl_resources [list "SLICE" "LUT" "FF" "DSP" "BRAM" "SRL"]
    set required_impl_resources_format "%8s %8s %8s %8s %8s %8s"

    set impl_timing_path "TimingReport"
    set required_impl_timing [list "TargetClockPeriod" "AchievedClockPeriod"]
    set required_impl_timing_format "%18s %18s"

    proc build_export_header { fd_out arch_param_names arch_param_names_format } {

        variable required_impl_resources
        variable required_impl_resources_format

        variable required_impl_timing
        variable required_impl_timing_format

        set line_format "${arch_param_names_format} ${required_impl_resources_format} ${required_impl_timing_format}"
        set table_header [concat ${arch_param_names} ${required_impl_resources} ${required_impl_timing}]
        
        puts [format $line_format {*}$table_header]
        puts $fd_out [format $line_format {*}$table_header]
    }

    proc build_export_line { fd_out fname arch_param arch_param_names arch_param_names_format} {
        variable impl_resources_path
        variable required_impl_resources
        variable required_impl_resources_format
        
        variable impl_timing_path
        variable required_impl_timing
        variable required_impl_timing_format

        set line_format "${arch_param_names_format} ${required_impl_resources_format} ${required_impl_timing_format}"

        if { [catch {open $fname} xmlfd] } {
            set nparams [expr [llength ${required_impl_resources}] + [llength ${required_impl_timing}]]
            for {set i 1} {$i <= $nparams} {incr i} {
                append results " FAILED"
            }
            set line [concat $arch_param $results]
            puts [format $line_format {*}$line]
            puts $fd_out [format $line_format {*}$line]
        } else {        
            set xmlcontent [read $xmlfd]
            close $xmlfd

            set xmldoc [dom parse $xmlcontent]
            set root [$xmldoc documentElement]

            set impl_resources_node [$root selectNodes ${impl_resources_path}]
            set impl_timing_node [$root selectNodes ${impl_timing_path}]

            set line [concat $arch_param                            \
                [get_fields ${required_impl_resources} ${impl_resources_node}]    \
                [get_fields ${required_impl_timing} ${impl_timing_node}]          ]
            puts [format $line_format {*}$line]
            puts $fd_out [format $line_format {*}$line]
        }
    }
}
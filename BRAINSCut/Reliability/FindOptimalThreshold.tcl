#!/usr/bin/tclsh
proc find_optimal {ANNResult ManaulTrace} {
    puts "proc"
      set ANNOutput [b2 load image $ANNResult data-type= float-single ]
      set Reference_image [b2 load image $ManaulTrace data-type= float-single ]
      set Reference [b2 threshold image $Reference_image 0.5 ]
    
    set optimavalue_for_RO 0
    set RO 0
    set optimavalue_for_SI 0
    set SI 0
    set i -0.05
    while { $i < 1.0 } {

      puts " /*"
      puts "  * iteration $i"
      puts "  */"

      set ANNCutOut [b2 threshold image $ANNOutput $i ]
      set AndMasks [b2 and masks $ANNCutOut $Reference ]
      set OrMasks [b2 or masks $ANNCutOut $Reference ]
      
      set Volume_AndMasks [lreplace [split [b2 measure volume mask $AndMasks] " "] 0 0 ]
      set Volume_AndMasks [lreplace [split $Volume_AndMasks "\\"] end end]
      puts "AND Volume = $Volume_AndMasks"
      set Volume_OrMasks  [lreplace [split [b2 measure volume mask $OrMasks] " "] 0 0 ]
      set Volume_OrMasks  [lreplace [split $Volume_OrMasks "\\"] end end]
      puts "OR Volume = $Volume_OrMasks"

      set Volume_Ref [lreplace [split [b2 measure volume mask $Reference] " "] 0 0 ]
      set Volume_ANN [lreplace [split [b2 measure volume mask $ANNCutOut] " "] 0 0 ]
      set Volume_Ref [lreplace [split $Volume_Ref  "\\"] end end]
      set Volume_ANN [lreplace [split $Volume_ANN  "\\"] end end]

      set temp_RO [expr $Volume_AndMasks / $Volume_OrMasks ]
      
      if {$temp_RO > $RO} {
        set RO $temp_RO
        set optimavalue_for_RO $i
      } 

      set A [expr 2 * $Volume_AndMasks]
      set B [expr $Volume_Ref + $Volume_ANN]
      set temp_SI [expr $A / $B ]

      if {$temp_SI > $SI} {
        set SI $temp_SI
        set optimavalue_for_SI $i
      } 
      puts " $i : $Volume_Ref : $Volume_ANN: $temp_RO : $temp_SI "
      set i [expr $i +0.01]
    
    }
    puts "Optimal value for RO is $RO at $optimavalue_for_RO"
    puts "Optimal value for SI is $SI at $optimavalue_for_SI"
  b2 destroy everything

}

set ann [lindex $argv 0]
set mannual [lindex $argv 1]

puts $ann
puts $mannual
if  { $argc < 2 }  { puts "USAGE:: [lindex $argv 0] ANN Manual" }
puts "start main"

puts "Th : Reference_volume : ANN_volume : RO : SI "
find_optimal  $ann $mannual





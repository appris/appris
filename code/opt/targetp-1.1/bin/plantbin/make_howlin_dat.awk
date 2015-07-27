# 
# FILE scripts/make_howlin_dat.awk
# usage 


/^/ {
                to_print=$0
                first_field=$1
                if (first_field=="l1no") { to_print=l1no }
                if (first_field=="l2no") { to_print=l2no }
                if (first_field=="l3no") { to_print=l3no }
                if (first_field=="sequencefile") { to_print=sequencefile }
                if (first_field=="synfile") { to_print=synfile number }   
                if (first_field=="noofseq") { to_print=noofseq }
                print to_print 
}

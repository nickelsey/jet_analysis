## install submit scripts


INSTALL( PROGRAMS qwrap.sh DESTINATION submit/ )
INSTALL( FILES pbs_submit_dijet_imbalance.py DESTINATION submit/ )
INSTALL( FILES pbs_submit_basic_dijet_imbalance.py DESTINATION submit/ )
INSTALL( FILES y6_bad_tower.txt DESTINATION submit/ )
INSTALL( FILES y7_y6_bad_tower.txt DESTINATION submit/ )
INSTALL( FILES y14_y6_bad_tower.txt DESTINATION submit/ )
INSTALL( FILES y14_bad_tower.txt DESTINATION submit/ )
INSTALL( FILES y14_bad_run.txt DESTINATION submit/ )
INSTALL( FILES y14_mb_file_list.txt DESTINATION submit/ )
INSTALL( DIRECTORY pp_list DESTINATION submit/
         FILES_MATCHING PATTERN "CMakeLists.txt" EXCLUDE
         PATTERN "._*" EXCLUDE
         PATTERN "*")

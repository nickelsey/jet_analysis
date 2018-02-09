## install submit scripts

install(PROGRAMS qwrap.sh DESTINATION submit)
install(FILES pbs_submit_dijet_imbalance.py DESTINATION submit)
install(FILES pbs_submit_basic_dijet_imbalance.py DESTINATION submit)
install(FILES pbs_submit_check_refmult.py DESTINATION submit)
install(FILES pbs_submit_comprehensive_dq.py DESTINATION submit)
install(FILES create_runid_list.py DESTINATION submit)
install(FILES y6_bad_tower.txt DESTINATION submit)
install(FILES y7_y6_bad_tower.txt DESTINATION submit)
install(FILES y14_y6_bad_tower.txt DESTINATION submit)
install(FILES y14_bad_tower.txt DESTINATION submit)
install(FILES y14_bad_run.txt DESTINATION submit)
install(FILES y14_mb_file_list.txt DESTINATION submit)
install(FILES y14_runids.root DESTINATION submit)
install(FILES y7_runids.root DESTINATION submit)
install(FILES empty_list.txt DESTINATION submit)
install(DIRECTORY pp_list DESTINATION submit
        FILES_MATCHING PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN "._*" EXCLUDE
        PATTERN "*")

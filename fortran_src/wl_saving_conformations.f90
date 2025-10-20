subroutine run_wl_low_energy( &
& n_flatness_checks, checks_flatness_every, &
& prng_seed, &
& n_monomers, sequence, interaction_matrix, &
& number_of_native_contacts, native_contacts, &
& n_energy_bins, energy_bin_edges, &
& flatness_trashold, &
& max_n_steps, &
& n_low_energy_structures, &
& n_conformations_per_bin, &
& save_conformation_every, &
& low_energies, &
& low_energies_contacts, &
& low_energies_n_contacts, &
& low_energies_n_native_contacts, &
& low_energy_bond_vectors, &
& bin_has_been_visited, &
& history_flat_histograms, &
& history_flat_ln_g_E, &
& history_f_factor, &
& saved_conformations, &
& saved_conformations_contacts, &
& saved_conformations_n_native_contacts)

    implicit none

    integer, parameter :: i8 = selected_int_kind(18)

    integer, intent(in) :: n_flatness_checks, checks_flatness_every
    integer, intent(in) :: n_monomers
    integer, intent(in):: sequence(n_monomers)
    double precision, intent(in) :: interaction_matrix(20,20)
    integer, intent(in) :: number_of_native_contacts
    integer, intent(in) :: native_contacts(number_of_native_contacts,2)
    integer, intent(in) :: prng_seed
    integer, intent(in) :: n_energy_bins
    integer, intent(in) :: n_low_energy_structures
    integer, intent(in) :: n_conformations_per_bin
    integer, intent(in) :: save_conformation_every
    double precision, intent(in) :: energy_bin_edges(n_energy_bins+1)
    double precision, intent(in) :: flatness_trashold
    integer(kind=i8), intent(in) :: max_n_steps
    
    integer, intent(out) :: history_flat_histograms(n_energy_bins,n_flatness_checks)
    double precision, intent(out) :: history_flat_ln_g_E(n_energy_bins,n_flatness_checks)
    logical, intent(out) :: bin_has_been_visited(n_energy_bins)
    double precision, intent(out) :: history_f_factor(n_flatness_checks)
    integer, intent(out) :: low_energies_n_contacts(n_low_energy_structures)
    integer, intent(out) :: low_energies_n_native_contacts(n_low_energy_structures)
    integer, intent(out) :: low_energies_contacts(n_low_energy_structures,2*n_monomers,2)
    integer, intent(out) :: low_energy_bond_vectors(n_low_energy_structures,n_monomers-1)
    double precision, intent(out) :: low_energies(n_low_energy_structures)
    integer, intent(out) :: saved_conformations(n_monomers-1,n_conformations_per_bin,n_energy_bins)
    integer, intent(out) :: saved_conformations_contacts(2*n_monomers,2,n_conformations_per_bin,n_energy_bins)
    integer, intent(out) :: saved_conformations_n_native_contacts(n_conformations_per_bin,n_energy_bins)
    double precision :: current_ln_g_E(n_energy_bins)

    integer(kind=i8) :: tot_n_steps
    integer :: tmp_obtained_bond_vectors(n_monomers-1)
    integer :: rot_tmp_obtained_bond_vectors(n_monomers-1)

    integer :: lattice_size
    double precision :: current_energy, test_energy 
    

    integer, allocatable :: current_conformation(:),test_conformation(:)

    double precision :: delta_energy,now_min_of_current_ln_g_E

    integer :: n_sites, i1, i2, i3, now_ind
    integer, allocatable :: ind_matrix(:,:,:), sub1_matrix(:), sub2_matrix(:), sub3_matrix(:)

    integer, dimension(6,3) :: lattice_versors
    integer, dimension(4,3) :: lattice_XY_versors,lattice_YZ_versors,lattice_XZ_versors
    integer, dimension(3) :: chosen_versor_crank
    integer :: neigh_along_versor
    integer, dimension(6) :: neigh_list

    integer :: i,i_low_energy


    integer :: i_step
    integer :: chosen_monomer, distance2, dx, dy, dz
    integer :: distance2_crank, dx_crank, dy_crank, dz_crank, rand_dir_crank
    logical :: there_is_overlap, is_accepted
    integer :: randi1

    integer :: h_moves_pro, h_moves_fea, h_moves_acc
    integer :: t_moves_pro, t_moves_fea, t_moves_acc
    integer :: c_moves_pro, c_moves_fea, c_moves_acc
    integer :: k_moves_pro, k_moves_fea, k_moves_acc

    integer :: current_number_of_contacts
    integer :: number_of_current_native_contacts
    integer, allocatable :: current_contacts(:,:)
    integer, allocatable :: current_native_contacts(:,:)

    integer :: now_bin

    double precision, allocatable :: energy_bin_centers(:)
    
    integer, allocatable :: energy_bin_heights(:)
    double precision :: now_f_modify

    integer :: n_seed
    integer, allocatable :: seed(:)

    logical :: is_duplicate
    integer :: i_low_energy_struct, i_bond,i_history_histograms

    integer :: current_n_bins
    double precision :: current_mean_hist_value,current_min_hist_value
    integer :: current_n_positive_flatness_checks

    integer :: i_step_save_conformations, rand_conf_ind
    integer :: n_saved_conformations_each_bin(n_energy_bins)

    !PRNG initialization
    call random_seed() 
    call random_seed(size=n_seed)
    allocate(seed(n_seed))
    seed = prng_seed
    call random_seed(put=seed)
    deallocate(seed)


    allocate(energy_bin_centers(n_energy_bins))
    allocate(energy_bin_heights(n_energy_bins))
    allocate(current_conformation(n_monomers))
    allocate(test_conformation(n_monomers))
    

    do i = 1, n_energy_bins
        energy_bin_centers(i) = 0.5*(energy_bin_edges(i) + energy_bin_edges(i+1))
    end do

    i_history_histograms=0
    current_energy = 0.0
    test_energy = 0.0
    delta_energy = 0.0
    now_f_modify = exp(1.0)
    current_ln_g_E = 0.0
    energy_bin_heights = 0
    bin_has_been_visited = .false.
    history_flat_histograms = 0
    history_flat_ln_g_E = 0.0
    history_f_factor = 0.0
    low_energies = 1.0e6
    low_energies_contacts = 0
    low_energies_n_contacts = 0
    low_energies_n_native_contacts = 0
    low_energy_bond_vectors = 0
    current_number_of_contacts = 0
    number_of_current_native_contacts = 0
    current_n_positive_flatness_checks = 0
    current_n_bins = 0
    current_mean_hist_value = 0.0
    current_min_hist_value = 0.0
    i_step = 0
    now_bin = 0
    saved_conformations = 0
    n_saved_conformations_each_bin = 0
    saved_conformations_contacts = -1
    saved_conformations_n_native_contacts = -1

    tot_n_steps = 0_i8

    neigh_list = 0

    h_moves_pro = 0
    h_moves_fea = 0
    h_moves_acc = 0
    t_moves_pro = 0
    t_moves_fea = 0
    t_moves_acc = 0
    c_moves_pro = 0
    c_moves_fea = 0
    c_moves_acc = 0
    k_moves_pro = 0
    k_moves_fea = 0
    k_moves_acc = 0


    lattice_size=n_monomers+10
    do i=1,n_monomers
        current_conformation(i) = 30*lattice_size+i+5
    end do
    test_conformation=current_conformation

    n_sites = lattice_size**3

    lattice_versors(1,:) = [-1, 0, 0]
    lattice_versors(2,:) = [1, 0, 0]
    lattice_versors(3,:) = [0, -1, 0]
    lattice_versors(4,:) = [0, 1, 0]
    lattice_versors(5,:) = [0, 0, -1]
    lattice_versors(6,:) = [0, 0, 1]

    lattice_XY_versors(1,:) = [-1, 0, 0]
    lattice_XY_versors(2,:) = [1, 0, 0]
    lattice_XY_versors(3,:) = [0, -1, 0]
    lattice_XY_versors(4,:) = [0, 1, 0]

    lattice_YZ_versors(1,:) = [0, -1, 0]
    lattice_YZ_versors(2,:) = [0, 1, 0]
    lattice_YZ_versors(3,:) = [0, 0, -1]
    lattice_YZ_versors(4,:) = [0, 0, 1]

    lattice_XZ_versors(1,:) = [-1, 0, 0]
    lattice_XZ_versors(2,:) = [1, 0, 0]
    lattice_XZ_versors(3,:) = [0, 0, -1]
    lattice_XZ_versors(4,:) = [0, 0, 1]

    allocate(ind_matrix(lattice_size, lattice_size, lattice_size))
    allocate(sub1_matrix(n_sites))
    allocate(sub2_matrix(n_sites))
    allocate(sub3_matrix(n_sites))
    now_ind=0
    do i1 = 1, lattice_size
        do i2 = 1, lattice_size
            do i3 = 1, lattice_size
                now_ind = now_ind + 1
                ind_matrix(i1, i2, i3) = now_ind
                sub1_matrix(now_ind) = i1
                sub2_matrix(now_ind) = i2
                sub3_matrix(now_ind) = i3
            end do
        end do
    end do
    allocate(current_contacts(2*n_monomers,2))
    allocate(current_native_contacts(number_of_native_contacts,2))
    current_contacts = 0
    current_native_contacts = 0


    call calculate_energy_of_conformation(lattice_size, &
    & ind_matrix,lattice_versors, &
    & n_monomers, current_conformation, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & sequence, interaction_matrix, current_energy)
    test_conformation=current_conformation
    test_energy=current_energy

    call return_current_contacts(lattice_size, &
    & ind_matrix,lattice_versors, &
    & n_monomers, current_conformation, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & current_number_of_contacts, current_contacts)

    call return_number_of_current_native_contacts(n_monomers,&
    & current_contacts, native_contacts, &
    & current_number_of_contacts, number_of_native_contacts, &
    & current_native_contacts, number_of_current_native_contacts)


    i_step = 0
    i_step_save_conformations = 0
    do 
        tot_n_steps = tot_n_steps + 1_i8
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!                                              !!
        !!            CHOOSE MOVE                       !!
        !!                                              !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! initialize test_conformation
        test_conformation=current_conformation

        ! choose a monomer
        call generate_randi(1,n_monomers,chosen_monomer)

        ! if chosen_monomer is 1, do a head move:
        if (chosen_monomer==1) then
            h_moves_pro = h_moves_pro + 1
            call return_neigh_list(lattice_size, &
            & sub1_matrix(current_conformation(2)), &
            & sub2_matrix(current_conformation(2)), &
            & sub3_matrix(current_conformation(2)), &
            & ind_matrix,lattice_versors,neigh_list)
            call generate_randi(1,6,randi1)
            test_conformation(1)=neigh_list(randi1)

            call check_for_overlaps_one_monomer(n_monomers,test_conformation, &
            & chosen_monomer, there_is_overlap)


            if (.not. there_is_overlap) then
                h_moves_fea = h_moves_fea + 1
                call mc_metropolis_one_monomer(lattice_size, &
                & ind_matrix, lattice_versors, &
                & n_monomers, current_conformation, test_conformation, &
                & chosen_monomer, &
                & sub1_matrix, sub2_matrix, sub3_matrix, &
                & sequence, interaction_matrix, delta_energy)
                test_energy = current_energy + delta_energy
                call compare_G_E(current_energy,test_energy,n_energy_bins,&
                & energy_bin_edges,current_ln_g_E,is_accepted)
                if (is_accepted) then
                    current_conformation=test_conformation
                    current_energy=current_energy+delta_energy
                    h_moves_acc = h_moves_acc + 1
                end if
            end if

            ! if chosen_monomer is n_monomers, do a tail move:
        else if (chosen_monomer==n_monomers) then
            t_moves_pro = t_moves_pro + 1
            call return_neigh_list(lattice_size, &
            & sub1_matrix(current_conformation(n_monomers-1)), &
            & sub2_matrix(current_conformation(n_monomers-1)), &
            & sub3_matrix(current_conformation(n_monomers-1)), &
            & ind_matrix,lattice_versors,neigh_list)
            call generate_randi(1,6,randi1)
            test_conformation(n_monomers)=neigh_list(randi1)

            call check_for_overlaps_one_monomer(n_monomers,test_conformation, &
            & chosen_monomer, there_is_overlap)

            if (.not. there_is_overlap) then
                t_moves_fea = t_moves_fea + 1
                call mc_metropolis_one_monomer(lattice_size, &
                & ind_matrix, lattice_versors, &
                & n_monomers, current_conformation, test_conformation, &
                & chosen_monomer, &
                & sub1_matrix, sub2_matrix, sub3_matrix, &
                & sequence, interaction_matrix, delta_energy)
                test_energy = current_energy + delta_energy
                call compare_G_E(current_energy,test_energy,n_energy_bins,&
                & energy_bin_edges,current_ln_g_E,is_accepted)
                if (is_accepted) then
                    current_conformation=test_conformation
                    current_energy=current_energy+delta_energy
                    t_moves_acc = t_moves_acc + 1
                end if
            end if

        else ! attempt corner flip

            c_moves_pro = c_moves_pro + 1

            ! check if it is feasible
            call calculate_distance2_pbc(lattice_size, &
            & current_conformation(chosen_monomer-1),current_conformation(chosen_monomer+1), &
            & sub1_matrix,sub2_matrix,sub3_matrix,distance2,dx,dy,dz)


            if (distance2 == 2) then ! corner flip is feasible

                ! check the plane in which to perform corner flip
                if (abs(dx) == 1 .and. abs(dz) == 1) then ! we are in the XZ plane
                    ! distinguish between two possible cases
                    if (sub1_matrix(current_conformation(chosen_monomer+1)) &
                    & == sub1_matrix(current_conformation(chosen_monomer))) then
                        test_conformation(chosen_monomer) = ind_matrix(sub1_matrix(current_conformation(chosen_monomer-1)),&
                        & sub2_matrix(current_conformation(chosen_monomer)), &
                        & sub3_matrix(current_conformation(chosen_monomer+1)))
                    else
                        test_conformation(chosen_monomer) = ind_matrix(sub1_matrix(current_conformation(chosen_monomer+1)),&
                        & sub2_matrix(current_conformation(chosen_monomer)),&
                        & sub3_matrix(current_conformation(chosen_monomer-1)))
                    end if
                else if (abs(dx) == 1 .and. abs(dy) == 1) then ! we are in the XY plane
                    if (sub1_matrix(current_conformation(chosen_monomer+1)) &
                    & == sub1_matrix(current_conformation(chosen_monomer))) then
                        test_conformation(chosen_monomer) = ind_matrix(sub1_matrix(current_conformation(chosen_monomer-1)),&
                        & sub2_matrix(current_conformation(chosen_monomer+1)), &
                        & sub3_matrix(current_conformation(chosen_monomer)))
                    else
                        test_conformation(chosen_monomer) = ind_matrix(sub1_matrix(current_conformation(chosen_monomer+1)),&
                        & sub2_matrix(current_conformation(chosen_monomer-1)),&
                        & sub3_matrix(current_conformation(chosen_monomer)))
                    end if
                else if (abs(dy) == 1 .and. abs(dz) == 1) then ! we are in the YZ plane
                    if (sub2_matrix(current_conformation(chosen_monomer+1)) &
                    & == sub2_matrix(current_conformation(chosen_monomer))) then
                        test_conformation(chosen_monomer) = ind_matrix(sub1_matrix(current_conformation(chosen_monomer)),&
                        & sub2_matrix(current_conformation(chosen_monomer-1)), &
                        & sub3_matrix(current_conformation(chosen_monomer+1)))
                    else
                        test_conformation(chosen_monomer) = ind_matrix(sub1_matrix(current_conformation(chosen_monomer)),&
                        & sub2_matrix(current_conformation(chosen_monomer+1)),&
                        & sub3_matrix(current_conformation(chosen_monomer-1)))
                    end if ! exhausted cases 1 or 2
                end if ! exhausted cases for different planes

                ! check for overlaps in the test conformation after corner flip
                call check_for_overlaps_one_monomer(n_monomers,test_conformation, &
                & chosen_monomer, there_is_overlap)

                if (.not. there_is_overlap) then
                    c_moves_fea = c_moves_fea + 1
                    call mc_metropolis_one_monomer(lattice_size, &
                    & ind_matrix, lattice_versors, &
                    & n_monomers, current_conformation, test_conformation, &
                    & chosen_monomer, &
                    & sub1_matrix, sub2_matrix, sub3_matrix, &
                    & sequence, interaction_matrix, delta_energy)
                    test_energy = current_energy + delta_energy
                    call compare_G_E(current_energy,test_energy,n_energy_bins,&
                    & energy_bin_edges,current_ln_g_E,is_accepted)
                    if (is_accepted) then
                        current_conformation=test_conformation
                        current_energy=current_energy+delta_energy
                        c_moves_acc = c_moves_acc + 1
                    end if
                end if

                if (there_is_overlap .and. (chosen_monomer < n_monomers-1)) then

                    k_moves_pro = k_moves_pro + 1

                    ! re-initialize test_conformation for crankshaft move
                    test_conformation=current_conformation

                    call calculate_distance2_pbc(lattice_size, &
                    & current_conformation(chosen_monomer - 1),current_conformation(chosen_monomer+2), &
                    & sub1_matrix,sub2_matrix,sub3_matrix,distance2_crank,dx_crank,dy_crank,dz_crank)

                    if (distance2_crank==1.0) then

                        call generate_randi(1,4,rand_dir_crank)

                        if (abs(dx_crank)==1) then
                            chosen_versor_crank = lattice_YZ_versors(rand_dir_crank,1:3)
                        else if (abs(dy_crank)==1) then
                            chosen_versor_crank = lattice_XZ_versors(rand_dir_crank,1:3)
                        else
                            chosen_versor_crank = lattice_XY_versors(rand_dir_crank,1:3)
                        end if

                        ! produce test conformation
                        call return_neigh_along_versor(lattice_size,&
                        & sub1_matrix(current_conformation(chosen_monomer-1)),&
                        & sub2_matrix(current_conformation(chosen_monomer-1)), &
                        & sub3_matrix(current_conformation(chosen_monomer-1)),&
                        & ind_matrix,chosen_versor_crank,neigh_along_versor)
                        test_conformation(chosen_monomer) = neigh_along_versor
                        call return_neigh_along_versor(lattice_size,&
                        & sub1_matrix(current_conformation(chosen_monomer+2)),&
                        & sub2_matrix(current_conformation(chosen_monomer+2)), &
                        & sub3_matrix(current_conformation(chosen_monomer+2)),&
                        & ind_matrix,chosen_versor_crank,neigh_along_versor)
                        test_conformation(chosen_monomer+1) = neigh_along_versor

                        ! call check_for_overlaps(n_monomers,test_conformation,there_is_overlap)
                        call check_for_overlaps_two_monomers(n_monomers,test_conformation, &
                        & chosen_monomer, there_is_overlap)


                        if (.not. there_is_overlap) then
                            k_moves_fea = k_moves_fea + 1
                            call mc_metropolis_two_monomers(lattice_size, &
                            & ind_matrix, lattice_versors, &
                            & n_monomers, current_conformation, test_conformation, &
                            & chosen_monomer, &
                            & sub1_matrix, sub2_matrix, sub3_matrix, &
                            & sequence, interaction_matrix, delta_energy)
                            test_energy = current_energy + delta_energy
                            call compare_G_E(current_energy,test_energy,n_energy_bins,&
                            & energy_bin_edges,current_ln_g_E,is_accepted)
                            if (is_accepted) then
                                current_conformation=test_conformation
                                current_energy=current_energy+delta_energy
                                k_moves_acc = k_moves_acc + 1
                            end if
                        end if

                    end if
                end if
            end if  ! close if (distance2 == 2), i.e. check for feasible cornerflip (which is also)
            ! a necessary condition for crank feasibility
        end if



        do now_bin = 1,n_energy_bins
            if (current_energy < energy_bin_edges(now_bin+1)) then
                exit
            end if
        end do
        
        energy_bin_heights(now_bin)=energy_bin_heights(now_bin)+1

        if (.not.(bin_has_been_visited(now_bin))) then
            
            now_min_of_current_ln_g_E=100000000000.0d0
            do i = 1, n_energy_bins
                if (bin_has_been_visited(i)) then
                    now_min_of_current_ln_g_E = min(&
                    & now_min_of_current_ln_g_E, &
                    & 1.0)
                end if
            end do
            
            if (now_min_of_current_ln_g_E<100000000000.0d0) then
                current_ln_g_E(now_bin)=now_min_of_current_ln_g_E
            else
                current_ln_g_E(now_bin)=current_ln_g_E(now_bin)+log(now_f_modify)
            end if
            bin_has_been_visited(now_bin)=.true.
            energy_bin_heights=0
        else
            current_ln_g_E(now_bin)=current_ln_g_E(now_bin)+log(now_f_modify)
        end if


        


        if (current_energy < low_energies(n_low_energy_structures)) then
            i_low_energy=n_low_energy_structures
            do i_low_energy=n_low_energy_structures-1,1,-1
                if (current_energy >= low_energies(i_low_energy)) then
                    exit
                end if
            end do
            i_low_energy=i_low_energy+1

            do i1=1,n_monomers-1
                call calculate_distance2_pbc(lattice_size, current_conformation(i1), current_conformation(i1+1), &
                    sub1_matrix, sub2_matrix, sub3_matrix, distance2, dx, dy, dz)
                if (dx==-1) then
                    tmp_obtained_bond_vectors(i1)=1
                else if (dx==1) then
                    tmp_obtained_bond_vectors(i1)=2
                else if (dy==-1) then
                    tmp_obtained_bond_vectors(i1)=3
                else if (dy==1) then
                    tmp_obtained_bond_vectors(i1)=4
                else if (dz==-1) then
                    tmp_obtained_bond_vectors(i1)=5
                else
                    tmp_obtained_bond_vectors(i1)=6
                end if

            end do
            call orient_chain(n_monomers-1, tmp_obtained_bond_vectors, rot_tmp_obtained_bond_vectors)





            is_duplicate = .false.
            do i_low_energy_struct = n_low_energy_structures,1,-1
                do i_bond = 1, n_monomers - 1
                    if (rot_tmp_obtained_bond_vectors(i_bond) /= low_energy_bond_vectors(i_low_energy_struct, i_bond)) exit
                end do
                if (i_bond > n_monomers - 1) then
                    is_duplicate = .true.
                    exit
                end if
            end do

            if (.not. is_duplicate) then
                do i1 = n_low_energy_structures - 1, i_low_energy, -1
                    low_energies(i1 + 1) = low_energies(i1)
                    low_energy_bond_vectors(i1 + 1, :) = low_energy_bond_vectors(i1, :)
                    low_energies_contacts(i1 + 1,:,:) = low_energies_contacts(i1,:,:)
                    low_energies_n_contacts(i1 + 1)=low_energies_n_contacts(i1)
                    low_energies_n_native_contacts(i1 + 1)=low_energies_n_native_contacts(i1)
                end do


                call return_current_contacts(lattice_size, &
                & ind_matrix,lattice_versors, &
                & n_monomers, current_conformation, &
                & sub1_matrix, sub2_matrix, sub3_matrix, &
                & current_number_of_contacts, current_contacts)

                call return_number_of_current_native_contacts(n_monomers,&
                & current_contacts, native_contacts, &
                & current_number_of_contacts, number_of_native_contacts, &
                & current_native_contacts, number_of_current_native_contacts)

                low_energies(i_low_energy) = current_energy
                low_energy_bond_vectors(i_low_energy,:) = rot_tmp_obtained_bond_vectors
                low_energies_contacts(i_low_energy,:,:)=current_contacts
                low_energies_n_contacts(i_low_energy)=current_number_of_contacts
                low_energies_n_native_contacts(i_low_energy)=number_of_current_native_contacts
            end if


        end if
       


        if (i_step > 0 .and. mod(i_step, checks_flatness_every) == 0) then
            i_step = 0
            current_mean_hist_value = 0.0
            current_n_bins = 0
            do i = 1, n_energy_bins
                if (bin_has_been_visited(i)) then
                    current_n_bins=current_n_bins+1
                    current_mean_hist_value=current_mean_hist_value+energy_bin_heights(i)
                end if
            end do
            current_mean_hist_value = current_mean_hist_value / dble(current_n_bins)

            current_min_hist_value = 100000000.0d0
            do i = 1, n_energy_bins
                if (bin_has_been_visited(i)) then
                    current_min_hist_value=min(&
                        & current_min_hist_value, &
                        & dble(energy_bin_heights(i)))
                end if
            end do
            
            if (current_min_hist_value > flatness_trashold*current_mean_hist_value) then
                current_n_positive_flatness_checks = current_n_positive_flatness_checks + 1
                history_f_factor(current_n_positive_flatness_checks)=now_f_modify
                history_flat_histograms(:,current_n_positive_flatness_checks)=energy_bin_heights
                history_flat_ln_g_E(:,current_n_positive_flatness_checks)=current_ln_g_E
                energy_bin_heights = 0
                now_f_modify=now_f_modify**0.5
            end if
        end if

        if (i_step_save_conformations > 0 .and. mod(i_step_save_conformations, save_conformation_every) == 0) then
            i_step_save_conformations = 0
        
            do i1=1,n_monomers-1
                call calculate_distance2_pbc(lattice_size, current_conformation(i1), current_conformation(i1+1), &
                    sub1_matrix, sub2_matrix, sub3_matrix, distance2, dx, dy, dz)
                if (dx==-1) then
                    tmp_obtained_bond_vectors(i1)=1
                else if (dx==1) then
                    tmp_obtained_bond_vectors(i1)=2
                else if (dy==-1) then
                    tmp_obtained_bond_vectors(i1)=3
                else if (dy==1) then
                    tmp_obtained_bond_vectors(i1)=4
                else if (dz==-1) then
                    tmp_obtained_bond_vectors(i1)=5
                else
                    tmp_obtained_bond_vectors(i1)=6
                end if

            end do
            call orient_chain(n_monomers-1, tmp_obtained_bond_vectors, rot_tmp_obtained_bond_vectors)

            call return_current_contacts(lattice_size, &
                & ind_matrix,lattice_versors, &
                & n_monomers, current_conformation, &
                & sub1_matrix, sub2_matrix, sub3_matrix, &
                & current_number_of_contacts, current_contacts)

            call return_number_of_current_native_contacts(n_monomers,&
                & current_contacts, native_contacts, &
                & current_number_of_contacts, number_of_native_contacts, &
                & current_native_contacts, number_of_current_native_contacts)

            if (n_saved_conformations_each_bin(now_bin)<n_conformations_per_bin) then
                n_saved_conformations_each_bin(now_bin)=n_saved_conformations_each_bin(now_bin)+1
                saved_conformations(:,n_saved_conformations_each_bin(now_bin),now_bin)=rot_tmp_obtained_bond_vectors
                saved_conformations_contacts(:,:,n_saved_conformations_each_bin(now_bin),now_bin)=current_contacts
                saved_conformations_n_native_contacts(n_saved_conformations_each_bin(now_bin),now_bin)=&
                    & number_of_current_native_contacts
            else
                call generate_randi(1,n_conformations_per_bin,rand_conf_ind)
                saved_conformations(:,rand_conf_ind,now_bin)=rot_tmp_obtained_bond_vectors
                saved_conformations_contacts(:,:,rand_conf_ind,now_bin)=current_contacts
                saved_conformations_n_native_contacts(rand_conf_ind,now_bin)=number_of_current_native_contacts
            end if
        end if

        
        i_step = i_step + 1
        i_step_save_conformations = i_step_save_conformations + 1
        if (current_n_positive_flatness_checks >= n_flatness_checks) goto 104
        if (tot_n_steps >= max_n_steps) goto 104 !10^10 max steps

    end do


104 deallocate(energy_bin_centers, energy_bin_heights)
    deallocate(current_conformation, test_conformation)
    deallocate(ind_matrix, sub1_matrix, sub2_matrix, sub3_matrix)
    deallocate(current_contacts, current_native_contacts)

contains




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compare_G_E(current_energy,test_energy,n_energy_bins,&
    & energy_bin_edges,current_ln_g_E,is_accepted)
        implicit none
        double precision, intent(in) :: current_energy,test_energy
        integer, intent(in) :: n_energy_bins
        double precision, intent(in) :: energy_bin_edges(n_energy_bins+1)
        double precision, intent(in) :: current_ln_g_E(n_energy_bins)
        logical, intent(out) :: is_accepted
        integer :: now_bin=0 ,test_bin=0
        double precision :: now_log_g_E,test_log_g_E,rand_metropolis

        is_accepted = .false.


        do now_bin=1,n_energy_bins
            if (current_energy < energy_bin_edges(now_bin+1)) then
                exit
            end if
        end do

        do test_bin=1,n_energy_bins
            if (test_energy < energy_bin_edges(test_bin+1)) then
                exit
            end if
        end do

        now_log_g_E=current_ln_g_E(now_bin)
        test_log_g_E=current_ln_g_E(test_bin)

        if (now_log_g_E>=test_log_g_E) then
            is_accepted = .true.
        else
            call random_number(rand_metropolis)
            if (now_log_g_E-test_log_g_E > log(rand_metropolis)) then
                is_accepted = .true.
            end if
        end if

    end subroutine compare_G_E

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine return_neigh_list(lattice_size,i1,i2,i3,ind_matrix,lattice_versors,neigh_list)
        implicit none
        integer, intent(in) :: lattice_size, i1, i2, i3
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(out), dimension(6) :: neigh_list
        integer :: n1,n2,n3
        integer :: i

        do i = 1, 6
            n1 = mod(i1 + lattice_versors(i, 1) + lattice_size - 1, lattice_size) + 1
            n2 = mod(i2 + lattice_versors(i, 2) + lattice_size - 1, lattice_size) + 1
            n3 = mod(i3 + lattice_versors(i, 3) + lattice_size - 1, lattice_size) + 1
            neigh_list(i) = ind_matrix(n1, n2, n3)
        end do
    end subroutine return_neigh_list


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine return_neigh_along_versor(lattice_size,i1,i2,i3,ind_matrix,versor,neigh_along_versor)
        implicit none
        integer, intent(in) :: lattice_size, i1, i2, i3
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: versor(3)
        integer, intent(out) :: neigh_along_versor
        integer :: n1,n2,n3

        n1 = mod(i1 + versor(1) + lattice_size - 1, lattice_size) + 1
        n2 = mod(i2 + versor(2) + lattice_size - 1, lattice_size) + 1
        n3 = mod(i3 + versor(3) + lattice_size - 1, lattice_size) + 1
        neigh_along_versor = ind_matrix(n1, n2, n3)

    end subroutine return_neigh_along_versor


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculate_distance2_pbc(lattice_size, ind_1, ind_2, &
        sub1_matrix, sub2_matrix, sub3_matrix, distance2, dx, dy, dz)
        implicit none
        integer, intent(in) :: lattice_size, ind_1, ind_2
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(out) :: distance2, dx, dy, dz
        integer :: raw_dx, raw_dy, raw_dz
        integer :: half_lattice

        half_lattice = lattice_size / 2

        raw_dx = sub1_matrix(ind_2) - sub1_matrix(ind_1)
        raw_dy = sub2_matrix(ind_2) - sub2_matrix(ind_1)
        raw_dz = sub3_matrix(ind_2) - sub3_matrix(ind_1)

        if (raw_dx > half_lattice) then
            dx = raw_dx - lattice_size
        else if (raw_dx < -half_lattice) then
            dx = raw_dx + lattice_size
        else
            dx = raw_dx
        end if

        if (raw_dy > half_lattice) then
            dy = raw_dy - lattice_size
        else if (raw_dy < -half_lattice) then
            dy = raw_dy + lattice_size
        else
            dy = raw_dy
        end if

        if (raw_dz > half_lattice) then
            dz = raw_dz - lattice_size
        else if (raw_dz < -half_lattice) then
            dz = raw_dz + lattice_size
        else
            dz = raw_dz
        end if

        distance2 = dx**2 + dy**2 + dz**2

    end subroutine calculate_distance2_pbc


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_for_overlaps_one_monomer(n_monomers,test_conformation, &
    & chosen_monomer, there_is_overlap)
        implicit none
        integer, intent(in) :: n_monomers, chosen_monomer
        integer, intent(in) :: test_conformation(n_monomers)
        logical, intent(out) :: there_is_overlap
        integer :: i

        there_is_overlap = .false.
        if (chosen_monomer==1) then
            do i=2, n_monomers
                if (test_conformation(chosen_monomer)==test_conformation(i)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do

        else if (chosen_monomer==n_monomers) then
            do i=1, n_monomers-1
                if (test_conformation(chosen_monomer)==test_conformation(i)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do

        else
            do i=1, chosen_monomer-1
                if (test_conformation(chosen_monomer)==test_conformation(i)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do
            if (.not. there_is_overlap) then
                do i=chosen_monomer+1,n_monomers
                    if (test_conformation(chosen_monomer)==test_conformation(i)) then
                        there_is_overlap = .true.
                        exit
                    end if
                end do
            end if
        end if

    end subroutine check_for_overlaps_one_monomer


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_for_overlaps_two_monomers(n_monomers,test_conformation, &
    & chosen_monomer, there_is_overlap)
        implicit none
        integer, intent(in) :: n_monomers, chosen_monomer
        integer, intent(in) :: test_conformation(n_monomers)
        logical, intent(out) :: there_is_overlap
        integer :: i

        there_is_overlap = .false.

        do i=1, chosen_monomer-1
            if (test_conformation(chosen_monomer)==test_conformation(i)) then
                there_is_overlap = .true.
                exit
            end if
        end do

        if (.not. there_is_overlap) then
            do i=chosen_monomer+1,n_monomers
                if (test_conformation(chosen_monomer)==test_conformation(i)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do
        end if

        if (.not. there_is_overlap) then
            do i=1, (chosen_monomer+1)-1
                if (test_conformation((chosen_monomer+1))==test_conformation(i)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do
        end if

        if (.not. there_is_overlap) then
            do i=(chosen_monomer+1)+1,n_monomers
                if (test_conformation((chosen_monomer+1))==test_conformation(i)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do
        end if

    end subroutine check_for_overlaps_two_monomers

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mc_metropolis_one_monomer(lattice_size, &
    & ind_matrix, lattice_versors, &
    & n_monomers, current_conformation, test_conformation, &
    & chosen_monomer, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & sequence, interaction_matrix, delta_energy)

        implicit none
        integer, intent(in) :: lattice_size, n_monomers, chosen_monomer
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(in) :: current_conformation(n_monomers), test_conformation(n_monomers)
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(in) :: sequence(n_monomers)
        double precision, intent(in) :: interaction_matrix(20,20)

        double precision, intent(out) :: delta_energy

        double precision :: current_local_energy = 0.0, test_local_energy = 0.0
        integer :: i_neigh, i_check
        integer, dimension(6) :: neigh_list

        current_local_energy = 0.0
        call return_neigh_list(lattice_size, &
        & sub1_matrix(current_conformation(chosen_monomer)), &
        & sub2_matrix(current_conformation(chosen_monomer)), &
        & sub3_matrix(current_conformation(chosen_monomer)), &
        & ind_matrix,lattice_versors,neigh_list)
        do i_neigh=1,6
            do i_check = 1, n_monomers
                if (neigh_list(i_neigh) == current_conformation(i_check)) then
                    current_local_energy=current_local_energy+interaction_matrix(sequence(chosen_monomer),sequence(i_check))
                end if
            end do
        end do

        test_local_energy = 0.0
        call return_neigh_list(lattice_size, &
        & sub1_matrix(test_conformation(chosen_monomer)), &
        & sub2_matrix(test_conformation(chosen_monomer)), &
        & sub3_matrix(test_conformation(chosen_monomer)), &
        & ind_matrix,lattice_versors,neigh_list)
        do i_neigh=1,6
            do i_check = 1, n_monomers
                if (neigh_list(i_neigh) == test_conformation(i_check)) then
                    test_local_energy=test_local_energy+interaction_matrix(sequence(chosen_monomer),sequence(i_check))
                end if
            end do
        end do

        delta_energy = test_local_energy-current_local_energy

    end subroutine mc_metropolis_one_monomer


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mc_metropolis_two_monomers(lattice_size, &
    & ind_matrix, lattice_versors, &
    & n_monomers, current_conformation, test_conformation, &
    & chosen_monomer, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & sequence, interaction_matrix, delta_energy)

        implicit none
        integer, intent(in) :: lattice_size, n_monomers, chosen_monomer
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(in) :: current_conformation(n_monomers), test_conformation(n_monomers)
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(in) :: sequence(n_monomers)
        double precision, intent(in) :: interaction_matrix(20,20)

        double precision, intent(out) :: delta_energy

        double precision :: current_local_energy = 0.0, test_local_energy = 0.0
        integer :: i_neigh, i_check
        integer, dimension(6) :: neigh_list

        current_local_energy = 0.0
        call return_neigh_list(lattice_size, &
        & sub1_matrix(current_conformation(chosen_monomer)), &
        & sub2_matrix(current_conformation(chosen_monomer)), &
        & sub3_matrix(current_conformation(chosen_monomer)), &
        & ind_matrix,lattice_versors,neigh_list)
        do i_neigh=1,6
            do i_check = 1, n_monomers
                if (neigh_list(i_neigh) == current_conformation(i_check)) then
                    current_local_energy=current_local_energy+interaction_matrix(sequence(chosen_monomer),sequence(i_check))
                end if
            end do
        end do
        call return_neigh_list(lattice_size, &
        & sub1_matrix(current_conformation(chosen_monomer+1)), &
        & sub2_matrix(current_conformation(chosen_monomer+1)), &
        & sub3_matrix(current_conformation(chosen_monomer+1)), &
        & ind_matrix,lattice_versors,neigh_list)
        do i_neigh=1,6
            do i_check = 1, n_monomers
                if (neigh_list(i_neigh) == current_conformation(i_check)) then
                    current_local_energy=current_local_energy+interaction_matrix(sequence(chosen_monomer+1),sequence(i_check))
                end if
            end do
        end do

        test_local_energy = 0.0
        call return_neigh_list(lattice_size, &
        & sub1_matrix(test_conformation(chosen_monomer)), &
        & sub2_matrix(test_conformation(chosen_monomer)), &
        & sub3_matrix(test_conformation(chosen_monomer)), &
        & ind_matrix,lattice_versors,neigh_list)
        do i_neigh=1,6
            do i_check = 1, n_monomers
                if (neigh_list(i_neigh) == test_conformation(i_check)) then
                    test_local_energy=test_local_energy+interaction_matrix(sequence(chosen_monomer),sequence(i_check))
                end if
            end do
        end do
        call return_neigh_list(lattice_size, &
        & sub1_matrix(test_conformation(chosen_monomer+1)), &
        & sub2_matrix(test_conformation(chosen_monomer+1)), &
        & sub3_matrix(test_conformation(chosen_monomer+1)), &
        & ind_matrix,lattice_versors,neigh_list)
        do i_neigh=1,6
            do i_check = 1, n_monomers
                if (neigh_list(i_neigh) == test_conformation(i_check)) then
                    test_local_energy=test_local_energy+interaction_matrix(sequence(chosen_monomer+1),sequence(i_check))
                end if
            end do
        end do

        delta_energy = test_local_energy-current_local_energy

    end subroutine mc_metropolis_two_monomers


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine return_current_contacts(lattice_size, &
    & ind_matrix,lattice_versors, &
    & n_monomers, conformation, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & current_number_of_contacts, current_contacts)
        implicit none
        integer, intent(in) :: n_monomers,lattice_size
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(in) :: conformation(n_monomers)
        integer, intent(out) :: current_number_of_contacts
        integer, intent(out) :: current_contacts(2*n_monomers,2)
        integer :: i_bond, i_neigh, i_check
        integer :: i1,i2,i3,i4,group_start,group_end
        integer,dimension(2) :: tmp
        integer, dimension(6) :: neigh_list

        current_contacts = 0
        current_number_of_contacts = 0

        do i_bond=1,n_monomers-1
            call return_neigh_list(lattice_size, &
            & sub1_matrix(conformation(i_bond)), &
            & sub2_matrix(conformation(i_bond)), &
            & sub3_matrix(conformation(i_bond)), &
            & ind_matrix,lattice_versors,neigh_list)
            do i_neigh=1,6
                do i_check = i_bond+2, n_monomers
                    if (neigh_list(i_neigh) == conformation(i_check)) then
                        current_number_of_contacts = current_number_of_contacts + 1
                        current_contacts(current_number_of_contacts,1) = i_bond
                        current_contacts(current_number_of_contacts,2) = i_check
                    end if
                end do
            end do
        end do

        ! do i1 = 1, current_number_of_contacts-1
        !     do i2 = i1+1, current_number_of_contacts
        !         if (current_contacts(i1,1) > current_contacts(i2,1)) then
        !             tmp = current_contacts(i1,:)
        !             current_contacts(i1,:) = current_contacts(i2,:)
        !             current_contacts(i2,:) = tmp
        !         end if
        !     end do
        ! end do

        i1 = 1
        do while (i1 <= current_number_of_contacts)
            group_start = i1

            ! Find end of group with same first column value
            do i2 = i1 + 1, current_number_of_contacts
                if (current_contacts(i2,1) /= current_contacts(group_start,1)) exit
            end do
            group_end = i2 - 1

            ! Sort group by second column using bubble sort
            do i3 = group_start, group_end - 1
                do i4 = i3 + 1, group_end
                    if (current_contacts(i3,2) > current_contacts(i4,2)) then
                        tmp = current_contacts(i3,:)
                        current_contacts(i3,:) = current_contacts(i4,:)
                        current_contacts(i4,:) = tmp
                    end if
                end do
            end do

            i1 = group_end + 1
        end do

    end subroutine return_current_contacts


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine return_number_of_current_native_contacts(n_monomers,&
    & current_contacts, native_contacts, &
    & current_number_of_contacts, number_of_native_contacts, &
    & current_native_contacts, number_of_current_native_contacts)

        implicit none
        integer, intent(in) :: n_monomers
        integer, intent(in) :: current_number_of_contacts, number_of_native_contacts
        integer, intent(in) :: current_contacts(2*n_monomers,2), native_contacts(number_of_native_contacts,2)
        integer, intent(out) :: current_native_contacts(number_of_native_contacts,2)
        integer, intent(out) :: number_of_current_native_contacts
        integer :: i1,i2

        number_of_current_native_contacts = 0
        current_native_contacts = 0

        do i1=1, current_number_of_contacts
            do i2=1, number_of_native_contacts
                ! write(104,*) i1,i2, current_contacts(i1, 1),current_contacts(i1, 2),&
                !     &    native_contacts(i2, 1), native_contacts(i2, 2)
                if (all(current_contacts(i1, :) == native_contacts(i2, :))) then

                    number_of_current_native_contacts=number_of_current_native_contacts+1
                    current_native_contacts(number_of_current_native_contacts,:)=current_contacts(i1, :)
                end if
            end do
        end do

    end subroutine return_number_of_current_native_contacts

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine generate_randi(min,max,randi)
        implicit none
        integer, intent(in) :: min, max
        integer, intent(out) :: randi
        double precision :: rand_real

        call random_number(rand_real)
        randi = min + int(rand_real * (max - min + 1)) ! [[min, max]] inclusive

    end subroutine generate_randi


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                      previous subroutines
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculate_energy_of_conformation(lattice_size, &
    & ind_matrix, lattice_versors, &
    & n_monomers, conformation, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & sequence, interaction_matrix, energy)
        implicit none
        integer, intent(in) :: lattice_size, n_monomers
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(in) :: conformation(n_monomers)
        integer, intent(in) :: sequence(n_monomers)
        double precision, intent(in) :: interaction_matrix(20,20)
        double precision, intent(out) :: energy
        integer :: i_bond, i_neigh, i_check
        integer, dimension(6) :: neigh_list


        energy=0.0
        do i_bond=1,n_monomers-1
            call return_neigh_list(lattice_size, &
            & sub1_matrix(conformation(i_bond)), &
            & sub2_matrix(conformation(i_bond)), &
            & sub3_matrix(conformation(i_bond)), &
            & ind_matrix,lattice_versors,neigh_list)
            do i_neigh=1,6
                do i_check = i_bond+2, n_monomers
                    if (neigh_list(i_neigh) == conformation(i_check)) then
                        energy=energy+interaction_matrix(sequence(i_bond),sequence(i_check))
                    end if
                end do
            end do
        end do

    end subroutine calculate_energy_of_conformation

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compare_current_and_test_conformations_crank(n_monomers,&
    & current_conformation,test_conformation,&
    & chosen_monomer,conformations_are_different)
        implicit none
        integer, intent(in) :: n_monomers, chosen_monomer
        integer, intent(in) :: current_conformation(n_monomers), test_conformation(n_monomers)
        logical, intent(out) :: conformations_are_different
        integer :: i

        conformations_are_different = .false.
        do i=chosen_monomer-1, chosen_monomer+2
            if (.not. (current_conformation(i) == test_conformation(i)) ) then
                conformations_are_different = .true.
                exit
            end if
            if (conformations_are_different) then
                exit
            end if
        end do

    end subroutine compare_current_and_test_conformations_crank



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_for_overlaps(n_monomers,test_conformation,there_is_overlap)
        implicit none
        integer, intent(in) :: n_monomers
        integer, intent(in) :: test_conformation(n_monomers)
        logical, intent(out) :: there_is_overlap
        integer :: i,j

        there_is_overlap = .false.
        do i=1,n_monomers-1
            do j=i+1,n_monomers
                if (test_conformation(i)==test_conformation(j)) then
                    there_is_overlap = .true.
                    exit
                end if
            end do
            if (there_is_overlap) then
                exit
            end if
        end do

    end subroutine check_for_overlaps


    subroutine orient_chain(n, raw, out)
        implicit none
        integer, intent(in)  :: n
        integer, intent(in)  :: raw(n)
        integer, intent(out) :: out(n)
        integer :: i, k, l, b2, op

        ! 1) first-bond -> 1 mapping (rot1)
        ! integer, parameter :: rot1(6,6) = reshape([ &
        !     1,2,3,4,5,6,   2,1,3,4,5,6,   4,3,1,2,5,6,  &
        !     3,4,2,1,5,6,   6,5,3,4,1,2,   5,6,3,4,2,1   &
        !     ], shape=[6,6] )
        ! integer :: rot1(6,6)

        integer :: rot1(6,6)
        integer :: rx(4,6)

        rot1(1,:) = (/ 1,2,3,4,5,6 /)
        rot1(2,:) = (/ 2,1,3,4,5,6 /)
        rot1(3,:) = (/ 4,3,1,2,5,6 /)
        rot1(4,:) = (/ 3,4,2,1,5,6 /)
        rot1(5,:) = (/ 6,5,3,4,1,2 /)
        rot1(6,:) = (/ 5,6,3,4,2,1 /)

        rx(1,:) = (/ 1,2,3,4,5,6 /)
        rx(2,:) = (/ 1,2,5,6,4,3 /)
        rx(3,:) = (/ 1,2,6,5,3,4 /)
        rx(4,:) = (/ 1,2,4,3,6,5 /)


        ! 2) rotations about +x: id, +90, -90, 180
        ! integer, parameter :: rx(4,6) = reshape([ &
        !     1,2,3,4,5,6,   1,2,5,6,4,3,   1,2,6,5,3,4,   1,2,4,3,6,5 &
        !     ], shape=[4,6] )

        ! --- apply the first-bond->1 rotation to get a "temporary" out[]:
        do i = 1,n
            out(i) = rot1(raw(1),raw(i))
        end do

        ! --- find the first index k>1 with a "turn":
        k = 0
        do i = 2, n
            if ( raw(i) /= raw(1) ) then
                k = i
                exit
            end if
        end do
        if (k == 0) return    ! straight-line walker-nothing more to do

        ! --- see where that bond landed under rot1:
        b2 = out(k)

        ! --- pick the x-rotation to carry b2->3:
        select case (b2)
          case (3); op = 1 
          case (4); op = 4
          case (5); op = 3
          case (6); op = 2
          case default
            print *, 'Error: unexpected bond ', b2
            stop
        end select

        ! --- apply that x-rotation to the whole chain
        do i = 1,n
            out(i) = rx(op, out(i))
        end do

        ! --- now enforce the third-bond-out-of-plane type 5:
        ! find first l>k with bond not in {1,2,3,4}
        l = 0
        do i = k+1, n
            if (out(i) > 4) then
                l = i
                exit
            end if
        end do

        ! if we found such an l, and it's +z (6), reflect in xy plane -> swap 5-6
        if (l /= 0 .and. out(l) == 6) then
            do i = 1,n
                if      (out(i) == 5) then
                    out(i) = 6
                else if (out(i) == 6) then
                    out(i) = 5
                end if
            end do
        end if

    end subroutine orient_chain



end subroutine run_wl_low_energy




! subroutine reset_wl_low_energy()
!     implicit none
!     ! Deallocate arrays if allocated
!     if (allocated(energy_bin_centers)) deallocate(energy_bin_centers)
!     if (allocated(energy_bin_heights)) deallocate(energy_bin_heights)
!     if (allocated(current_conformation)) deallocate(current_conformation)
!     if (allocated(test_conformation)) deallocate(test_conformation)
!     if (allocated(ind_matrix)) deallocate(ind_matrix)
!     if (allocated(sub1_matrix)) deallocate(sub1_matrix)
!     if (allocated(sub2_matrix)) deallocate(sub2_matrix)
!     if (allocated(sub3_matrix)) deallocate(sub3_matrix)
!     if (allocated(current_contacts)) deallocate(current_contacts)
!     if (allocated(current_native_contacts)) deallocate(current_native_contacts)
!     ! Reset RNG to default state
!     call random_seed()
! end subroutine reset_wl_low_energy
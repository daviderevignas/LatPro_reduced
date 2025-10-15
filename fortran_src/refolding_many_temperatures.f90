subroutine run_refolding_simulations(n_refolding_steps, n_refolding_savings, n_monomers, &
    n_indepentend_refolding_runs_each_temperature, prng_seed, &
    interaction_matrix, number_of_native_contacts, native_contacts, sequence, &
    n_temperatures_to_do, temperatures_to_do, &
    n_initial_saw_steps, &
    all_rep_obtained_contacts, all_rep_obtained_energies, &
    all_rep_obtained_natives, all_rep_obtained_distances2, &
    all_rep_obtained_bond_vectors)
    implicit none

    ! Inputs
    integer, intent(in) :: n_refolding_steps, n_refolding_savings, n_monomers
    integer, intent(in) :: n_indepentend_refolding_runs_each_temperature, prng_seed
    double precision, intent(in) :: interaction_matrix(20, 20)
    integer, intent(in) :: number_of_native_contacts
    integer, intent(in) :: native_contacts(number_of_native_contacts, 2)
    integer, intent(in) :: sequence(n_monomers)
    integer, intent(in) :: n_temperatures_to_do
    double precision, intent(in) :: temperatures_to_do(n_temperatures_to_do)
    integer, intent(in) :: n_initial_saw_steps


    ! Outputs
    integer, intent(out) :: all_rep_obtained_natives(n_refolding_savings,n_indepentend_refolding_runs_each_temperature,n_temperatures_to_do)
    double precision, intent(out) :: all_rep_obtained_energies(n_refolding_savings,&
    & n_indepentend_refolding_runs_each_temperature,n_temperatures_to_do)
    integer, intent(out) :: all_rep_obtained_distances2(n_refolding_savings, &
    & (n_monomers*(n_monomers-1))/2, n_indepentend_refolding_runs_each_temperature,n_temperatures_to_do)
    integer, intent(out) :: all_rep_obtained_contacts(n_refolding_savings,2*n_monomers,&
    & 2,n_indepentend_refolding_runs_each_temperature,n_temperatures_to_do)
    integer, intent(out) :: all_rep_obtained_bond_vectors(n_refolding_savings,n_monomers-1,&
    & n_indepentend_refolding_runs_each_temperature,n_temperatures_to_do)

    ! INTEGER, ALLOCATABLE, INTENT(OUT) :: all_rep_obtained_natives(:,:)
    ! INTEGER, ALLOCATABLE, INTENT(OUT) :: all_rep_obtained_distances2(:,:,:)

    ! Internals
    integer :: now_obtained_contacts(n_refolding_savings,2*n_monomers,2)

    integer, allocatable :: now_n_obtained_natives(:)
    double precision, allocatable :: now_obtained_energies(:)
    integer, allocatable :: now_obtained_distances2(:,:)
    integer, allocatable :: now_obtained_bond_vectors(:,:)

    integer :: i_replica,i_T

    integer :: n_seed
    integer, allocatable :: seed(:)

    ! all_rep_obtained_contacts=-1
    ! all_rep_obtained_natives=-1
    ! all_rep_obtained_distances2=-1


    !PRNG initialization
    call random_seed(size=n_seed)
    allocate(seed(n_seed))
    seed = prng_seed
    call random_seed(put=seed)
    deallocate(seed)

    ! allocate(all_rep_obtained_natives(n_refolding_savings,n_indepentend_refolding_runs_each_temperature))
    ! allocate(all_rep_obtained_distances2(n_refolding_savings, (n_monomers*(n_monomers-1))/2, n_indepentend_refolding_runs_each_temperature))


    allocate(now_n_obtained_natives(n_refolding_savings))
    allocate(now_obtained_energies(n_refolding_savings))
    allocate(now_obtained_distances2(n_refolding_savings,(n_monomers*(n_monomers-1))/2))
    allocate(now_obtained_bond_vectors(n_refolding_savings,n_monomers-1))

    ! open(unit=104, file="debug_1.txt", status="unknown", action="write")
    do i_T=1,n_temperatures_to_do
        do i_replica=1,n_indepentend_refolding_runs_each_temperature
            call refolding(&
            & n_refolding_steps,n_refolding_savings,n_monomers,temperatures_to_do(i_T), &
            & sequence, interaction_matrix, &
            & number_of_native_contacts, native_contacts, &
            & n_initial_saw_steps, &
            & now_obtained_contacts, now_obtained_energies, &
            & now_n_obtained_natives,now_obtained_distances2, &
            & now_obtained_bond_vectors)

            all_rep_obtained_contacts(:,:,:,i_replica,i_T) = now_obtained_contacts
            all_rep_obtained_energies(:,i_replica,i_T) = now_obtained_energies
            all_rep_obtained_natives(:,i_replica,i_T) = now_n_obtained_natives
            all_rep_obtained_distances2(:,:,i_replica,i_T) = now_obtained_distances2
            all_rep_obtained_bond_vectors(:,:,i_replica,i_T) = now_obtained_bond_vectors
        end do
    end do

    ! close(104)
contains

    subroutine refolding( &
    & n_refolding_steps, n_refolding_savings,n_monomers,refolding_temperature, &
    & sequence, interaction_matrix, &
    & number_of_native_contacts, native_contacts, &
    & n_initial_saw_steps, &
    & obtained_contacts, obtained_energies, &
    & n_obtained_natives, obtained_distances2, &
    & obtained_bond_vectors)

        implicit none
        integer, intent(in) :: n_refolding_steps, n_refolding_savings, n_monomers
        double precision, intent(in) :: refolding_temperature
        integer, intent(in):: sequence(n_monomers)
        double precision, intent(in) :: interaction_matrix(20,20)

        integer, intent(in) :: number_of_native_contacts
        integer, intent(in) :: native_contacts(number_of_native_contacts,2)
        integer, intent(in) :: n_initial_saw_steps

        integer, intent(out) :: n_obtained_natives(n_refolding_savings)
        double precision, intent(out) :: obtained_energies(n_refolding_savings)
        integer, intent(out) :: obtained_distances2(n_refolding_savings, (n_monomers*(n_monomers-1))/2 )
        integer, intent(out) :: obtained_contacts(n_refolding_savings,2*n_monomers,2)

        ! integer :: current_bond_vectors(n_monomers-1)

        integer, intent(out) :: obtained_bond_vectors(n_refolding_savings,n_monomers-1)

        integer :: tmp_obtained_bond_vectors(n_monomers-1)


        integer :: lattice_size
        double precision :: current_energy=0.0, test_energy=0.0

        integer, allocatable :: current_conformation(:),test_conformation(:)

        double precision :: delta_energy

        integer :: n_sites, i1, i2, i3, now_ind
        integer, allocatable :: ind_matrix(:,:,:), sub1_matrix(:), sub2_matrix(:), sub3_matrix(:)

        integer, dimension(6,3) :: lattice_versors
        integer, dimension(4,3) :: lattice_XY_versors,lattice_YZ_versors,lattice_XZ_versors
        integer, dimension(3) :: chosen_versor_crank
        integer :: neigh_along_versor
        integer, dimension(6) :: neigh_list = 0

        integer :: i


        integer :: i_step, chosen_monomer, distance2, dx, dy, dz
        integer :: distance2_crank, dx_crank, dy_crank, dz_crank, rand_dir_crank
        logical :: there_is_overlap, is_accepted
        ! double precision :: now_manhattan_ete
        integer :: randi1

        integer :: h_moves_pro=0,h_moves_fea=0,h_moves_acc=0
        integer :: t_moves_pro=0,t_moves_fea=0,t_moves_acc=0
        integer :: c_moves_pro=0,c_moves_fea=0,c_moves_acc=0
        integer :: k_moves_pro=0,k_moves_fea=0,k_moves_acc=0

        integer :: current_number_of_contacts = 0
        integer :: number_of_current_native_contacts = 0
        integer, allocatable :: current_contacts(:,:)
        integer, allocatable :: current_native_contacts(:,:)


        integer :: save_every

        save_every=n_refolding_steps/n_refolding_savings


        allocate(current_conformation(n_monomers))
        allocate(test_conformation(n_monomers))


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


        

        do i_step=1,n_initial_saw_steps

            ! initialize test_conformation
            test_conformation=current_conformation

            ! choose a monomer
            call generate_randi(1,n_monomers,chosen_monomer)

            ! if chosen_monomer is 1, do a head move:
            if (chosen_monomer==1) then
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
                    current_conformation=test_conformation
                end if

                ! if chosen_monomer is n_monomers, do a tail move:
            else if (chosen_monomer==n_monomers) then
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
                    current_conformation=test_conformation
                end if

            else ! attempt corner flip

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

                    call check_for_overlaps_one_monomer(n_monomers,test_conformation, &
                    & chosen_monomer, there_is_overlap)
                    if (.not. there_is_overlap) then
                        current_conformation=test_conformation
                    end if

                    if (there_is_overlap .and. (chosen_monomer < n_monomers-1)) then
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
                            call check_for_overlaps_two_monomers(n_monomers,test_conformation, &
                            & chosen_monomer, there_is_overlap)
                            if (.not. there_is_overlap) then
                                current_conformation=test_conformation
                            end if
                        end if
                    end if
                end if  
            end if
        end do

        test_conformation=current_conformation

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



        do i_step = 0,n_refolding_steps-1

            ! SAVING:
            if (mod(i_step, save_every) == 0) then
                ! write(104,*) "================================================"
                ! write(104,*) i_step

                i3=0
                do i1=1,n_monomers-1
                    do i2=i1+1,n_monomers
                        i3=i3+1
                        call calculate_distance2_pbc(lattice_size, current_conformation(i1), current_conformation(i2), &
                            sub1_matrix, sub2_matrix, sub3_matrix, distance2, dx, dy, dz)
                        obtained_distances2(i_step/save_every+1,i3) = dx**2 + dy**2 + dz**2

                    end do
                end do

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
                call orient_chain(n_monomers-1, tmp_obtained_bond_vectors, obtained_bond_vectors(i_step/save_every+1,:))
                ! obtained_bond_vectors(i_step/save_every+1,:)=tmp_obtained_bond_vectors
                ! write(104,*) current_bond_vectors



                call return_current_contacts(lattice_size, &
                & ind_matrix,lattice_versors, &
                & n_monomers, current_conformation, &
                & sub1_matrix, sub2_matrix, sub3_matrix, &
                & current_number_of_contacts, current_contacts)

                ! do i1=1,current_number_of_contacts
                !      write(104,*) current_contacts(i1,:)
                ! end do




                call return_number_of_current_native_contacts(n_monomers,&
                & current_contacts, native_contacts, &
                & current_number_of_contacts, number_of_native_contacts, &
                & current_native_contacts, number_of_current_native_contacts)


                ! write(104,*) "----------"
                ! write(104,*) current_contacts
                ! write(104,*) native_contacts
                ! write(104,*) current_number_of_contacts
                ! write(104,*) number_of_native_contacts
                ! write(104,*) current_native_contacts
                ! write(104,*) number_of_current_native_contacts

                ! do i1=current_number_of_contacts+1,2*n_monomers
                !     current_contacts(i1,1)=-1
                !     current_contacts(i1,2)=-1
                ! end do
                ! write(104,*) current_energy

                obtained_contacts(i_step/save_every+1,:,:) = current_contacts
                obtained_energies(i_step/save_every+1) = current_energy
                n_obtained_natives(i_step/save_every+1) = number_of_current_native_contacts

                h_moves_pro=0
                h_moves_fea=0
                h_moves_acc=0

                t_moves_pro=0
                t_moves_fea=0
                t_moves_acc=0

                c_moves_pro=0
                c_moves_fea=0
                c_moves_acc=0

                k_moves_pro=0
                k_moves_fea=0
                k_moves_acc=0


            end if





            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!                                              !!
            !!            CHOOSE MOVE                       !!
            !!                                              !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! initialize test_conformation
            test_conformation=current_conformation

            ! choose a monomer
            call generate_randi(1,n_monomers,chosen_monomer)

            !write(debug_unit, *) chosen_monomer

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
                    & refolding_temperature, &
                    & sequence, interaction_matrix, is_accepted, delta_energy)
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
                    & refolding_temperature, &
                    & sequence, interaction_matrix, is_accepted, delta_energy)
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
                        & refolding_temperature, &
                        & sequence, interaction_matrix, is_accepted, delta_energy)
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
                                & refolding_temperature, &
                                & sequence, interaction_matrix, is_accepted, delta_energy)
                                !write(debug_unit, *) delta_energy
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




            ! now_manhattan_ete=0
            ! do i=1,n_monomers-1
            !     call calculate_distance2_pbc(lattice_size,current_conformation(i),current_conformation(i+1), &
            !         & sub1_matrix,sub2_matrix,sub3_matrix,distance2,dx, dy, dz)
            !     now_manhattan_ete=now_manhattan_ete + dble(abs(dx)+abs(dy)+abs(dz))
            ! end do
            ! if (now_manhattan_ete > dble(n_monomers-1)) then
            !     print *, 'PROBLEM WITH now_manhattan_ete!!'
            !     exit
            ! end if

        end do



    end subroutine refolding

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
    & refolding_temperature, &
    & sequence, interaction_matrix, is_accepted, delta_energy)

        implicit none
        integer, intent(in) :: lattice_size, n_monomers, chosen_monomer
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(in) :: current_conformation(n_monomers), test_conformation(n_monomers)
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(in) :: sequence(n_monomers)
        double precision, intent(in) :: interaction_matrix(20,20)
        double precision, intent(in) :: refolding_temperature

        logical, intent(out) :: is_accepted
        double precision, intent(out) :: delta_energy

        double precision :: current_local_energy = 0.0, test_local_energy = 0.0
        integer :: i_neigh, i_check
        integer, dimension(6) :: neigh_list
        double precision :: rand_metropolis

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
        is_accepted = .false.
        if (delta_energy<=0) then
            is_accepted = .true.
        else
            call random_number(rand_metropolis)
            if (exp(-delta_energy/refolding_temperature)>rand_metropolis) then
                is_accepted = .true.
            end if
        end if

    end subroutine mc_metropolis_one_monomer


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mc_metropolis_two_monomers(lattice_size, &
    & ind_matrix, lattice_versors, &
    & n_monomers, current_conformation, test_conformation, &
    & chosen_monomer, &
    & sub1_matrix, sub2_matrix, sub3_matrix, &
    & refolding_temperature, &
    & sequence, interaction_matrix, is_accepted, delta_energy)

        implicit none
        integer, intent(in) :: lattice_size, n_monomers, chosen_monomer
        integer, intent(in) :: ind_matrix(lattice_size,lattice_size,lattice_size)
        integer, intent(in) :: lattice_versors(6,3)
        integer, intent(in) :: current_conformation(n_monomers), test_conformation(n_monomers)
        integer, intent(in) :: sub1_matrix(lattice_size**3), sub2_matrix(lattice_size**3), sub3_matrix(lattice_size**3)
        integer, intent(in) :: sequence(n_monomers)
        double precision, intent(in) :: interaction_matrix(20,20)
        double precision, intent(in) :: refolding_temperature

        logical, intent(out) :: is_accepted
        double precision, intent(out) :: delta_energy

        double precision :: current_local_energy = 0.0, test_local_energy = 0.0
        integer :: i_neigh, i_check
        integer, dimension(6) :: neigh_list
        double precision :: rand_metropolis

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
        is_accepted = .false.
        if (delta_energy<=0) then
            is_accepted = .true.
        else
            call random_number(rand_metropolis)
            if (exp(-delta_energy/refolding_temperature)>rand_metropolis) then
                is_accepted = .true.
            end if
        end if

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
        integer :: i_mon, i_neigh, i_check
        integer :: i1,i2,i3,i4,group_start,group_end
        integer,dimension(2) :: tmp
        integer, dimension(6) :: neigh_list

        current_contacts = 0
        current_number_of_contacts = 0

        do i_mon=1,n_monomers-1
            call return_neigh_list(lattice_size, &
            & sub1_matrix(conformation(i_mon)), &
            & sub2_matrix(conformation(i_mon)), &
            & sub3_matrix(conformation(i_mon)), &
            & ind_matrix,lattice_versors,neigh_list)
            do i_neigh=1,6
                do i_check = i_mon+2, n_monomers
                    if (neigh_list(i_neigh) == conformation(i_check)) then
                        current_number_of_contacts = current_number_of_contacts + 1
                        current_contacts(current_number_of_contacts,1) = i_mon
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

            ! Find end of group with same first-column value
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
        integer :: i_mon, i_neigh, i_check
        integer, dimension(6) :: neigh_list


        energy=0.0
        do i_mon=1,n_monomers-1
            call return_neigh_list(lattice_size, &
            & sub1_matrix(conformation(i_mon)), &
            & sub2_matrix(conformation(i_mon)), &
            & sub3_matrix(conformation(i_mon)), &
            & ind_matrix,lattice_versors,neigh_list)
            do i_neigh=1,6
                do i_check = i_mon+2, n_monomers
                    if (neigh_list(i_neigh) == conformation(i_check)) then
                        energy=energy+interaction_matrix(sequence(i_mon),sequence(i_check))
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

        ! 1) first‐bond → 1 mapping (rot1)
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


        ! 2) rotations about +x: id, +90°, -90°, 180°
        ! integer, parameter :: rx(4,6) = reshape([ &
        !     1,2,3,4,5,6,   1,2,5,6,4,3,   1,2,6,5,3,4,   1,2,4,3,6,5 &
        !     ], shape=[4,6] )

        ! --- apply the first‐bond→1 rotation to get a “temporary” out[]:
        do i = 1,n
            out(i) = rot1(raw(1),raw(i))
        end do

        ! --- find the first index k>1 with a “turn”:
        k = 0
        do i = 2, n
            if ( raw(i) /= raw(1) ) then
                k = i
                exit
            end if
        end do
        if (k == 0) return    ! straight‐line walker—nothing more to do

        ! --- see where that bond landed under rot1:
        b2 = out(k)

        ! --- pick the x‐rotation to carry b2→3:
        select case (b2)
          case (3); op = 1     ! already −y
          case (4); op = 4     ! 180° about +x: +y→−y
          case (5); op = 3     ! −90° about +x: −z→−y
          case (6); op = 2     ! +90° about +x: +z→−y
          case default
            print *, 'Error: unexpected bond ', b2
            stop
        end select

        ! --- apply that x‐rotation to the whole chain
        do i = 1,n
            out(i) = rx(op, out(i))
        end do

        ! --- now enforce the third‐bond‐out‐of‐plane ⇒ type 5:
        ! find first l>k with bond ∉{1,2,3,4}
        l = 0
        do i = k+1, n
            if (out(i) > 4) then
                l = i
                exit
            end if
        end do

        ! if we found such an l, and it’s +z (6), reflect in xy plane → swap 5↔6
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


end subroutine run_refolding_simulations

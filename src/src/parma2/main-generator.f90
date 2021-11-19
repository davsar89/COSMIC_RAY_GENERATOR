
subroutine gen_parma_cr(seed, type, emin, emax, nb, glat, glong, alti, ener_out, u_out, v_out, w_out)
   use :: iso_c_binding
!  Generate cosmic-ray based on PARMA model
   implicit none
   integer, intent(in) :: nb
   integer, intent(in) :: seed
   real(8), intent(in) :: type
   real(8), intent(in) :: emin, emax
   real(8), intent(in) :: glat  ! Latitude (deg), -90 =< glat =< 90
   real(8), intent(in) :: glong ! Longitude (deg), -180 =< glong =< 180
   real(8), intent(in) :: alti ! altitude, km
   real(8), intent(out), dimension(nb) :: ener_out
   real(8), intent(out), dimension(nb) :: u_out
   real(8), intent(out), dimension(nb) :: v_out
   real(8), intent(out), dimension(nb) :: w_out
   real(8) :: u, v, w, e
   integer :: npart, nebin, nabin, i, ia, iday
   integer :: idummy, ie, imonth, ip, iyear, nevent
   integer:: IangPart
   real(8) :: amin, amax, ahigh, amid, astep, atable, etable, ehigh, emid, estep
   real(8) :: d, elog, g
   real(8) :: flux, totalflux
   real(8) :: gethp, getr, getd, getspecangfinal, getGeneration, getspec, grnd
   real(8) :: phi, r, s, sx
   real(8) :: PPII = acos(-1.0)
   parameter(npart=33) ! number of applicable particle
   parameter(nebin=1000) ! number of energy mesh (divided by log)
   parameter(nabin=100) ! number of angle mesh (divided by linear)
   dimension IangPart(0:npart)
   dimension ehigh(0:nebin), emid(nebin) ! higher and middle point of energy bin
   dimension ahigh(0:nabin), amid(nabin) ! higher and middle point of angular bin
   dimension etable(0:nebin)            ! probability table (0.0 for 0, 1.0 for nebin)
   dimension atable(0:nabin, 0:nebin)    ! probability table (0.0 for 0, 1.0 for nabin)
   dimension flux(0:nabin, nebin) ! Monte Carlo generated flux
   data IangPart/1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 6/ ! Angular particle ID

! Set number of particle to be generated, and initial random seed
   nevent = nb ! number of particles to be generated

   call sgrnd(seed) ! initial seed of random number (you can change the value)

! Set Particle ID
   ip = 0 ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
   if (IangPart(ip) .eq. 0) then
      write (*, *) 'Angular distribution is not available for the particle'
      stop
   end if

   if (emin .lt. 0.001) then
      write (*, *), "Energies below 1 MeV are not allowed."
      call abort
   end if

! Set Conditions (location, time, and local geometry)
   iyear = 2020  ! Year
   imonth = 9    ! Month
   iday = 10      ! Date
!glat=78.0   ! Latitude (deg), -90 =< glat =< 90
!glong=0.0 ! Longitude (deg), -180 =< glong =< 180
!alti=50.0   ! Altitude (km)
   g = 0.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin

   s = getHP(iyear, imonth, iday, idummy) ! Solar activity (W index) of the day
   r = getr(glat, glong)                ! Vertical cut-off rigidity (GV)
   d = getd(alti, glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

   if (isnan(r)) then
      write (*, *), "ERROR : rr is nan. Aborting."
      call abort
   end if
   if (isnan(s)) then
      write (*, *), "ERROR : s is nan. Aborting."
      call abort
   end if
   if (isnan(d)) then
      write (*, *), "ERROR : d is nan. Aborting."
      call abort
   end if

! Set energy and angle ranges for generation
!emin = 1.0e0  ! Minimum energy of particle, GeV
!emax = 1.0e3  ! Maximum energy of particle
   amin = -1.0   ! Minimum cosine of particle
   amax = 1.0   ! Maximum cosine of particle
!if(ip.eq.0.and.emin.lt.1.0e-8) emin=1.0e-8
!if(ip.ne.0.and.emin.lt.1.0e-2) emin=1.0e-2

! Make energy and angle mesh
   elog = log10(emin)
   estep = (log10(emax) - log10(emin))/nebin
   do ie = 0, nebin
      ehigh(ie) = 10d0**elog
      if (ie .ne. 0) emid(ie) = sqrt(ehigh(ie)*ehigh(ie - 1))
      elog = elog + estep
   end do

   astep = (amax - amin)/nabin
   do ia = 0, nabin
      ahigh(ia) = amin + astep*ia
      if (ia .ne. 0) amid(ia) = (ahigh(ia) + ahigh(ia - 1))*0.5
   end do

! Make probability table (absolute value)
   atable(:, :) = 0.0d0 ! initialization
   etable(:) = 0.0d0
   do ie = 1, nebin
      do ia = 1, nabin
         atable(ia, ie) = atable(ia - 1, ie) + getSpec(ip, s, r, d, emid(ie), g) &
                         & *getSpecAngFinal(iangpart(ip), s, r, d, emid(ie), g, amid(ia))*(2.0*PPII)*(ahigh(ia) - ahigh(ia - 1)) ! angular integrated value
         if (isnan(atable(ia, ie))) then
            write (*, *) "atable(ia,ie) is nan. Aborting."
            call abort
         end if
      end do
   end do
   do ie = 1, nebin
      etable(ie) = etable(ie - 1) + atable(nabin, ie)*(ehigh(ie) - ehigh(ie - 1)) ! energy integrated value
      if (isnan(etable(ie))) then
         write (*, *) "etable(ie) is nan. Aborting."
         call abort
      end if
   end do
   TotalFlux = etable(nebin) ! Total Flux (/cm2/s), used for normalization

! Make probability table (normalized to 1)
   do ie = 1, nebin
      etable(ie) = etable(ie)/etable(nebin)
      do ia = 1, nabin
         atable(ia, ie) = atable(ia, ie)/atable(nabin, ie)
         if (isnan(atable(ia, ie))) then
            write (*, *) "atable(ia,ie) is nan. Aborting."
            call abort
         end if
      end do
   end do

! Particle Generation

   do i = 1, nevent
      e = getGeneration(ie, nebin, ehigh(0), etable(0))    ! energy
      phi = 2.0*PPII*(grnd() - 0.5) ! azimuth angle (rad)
      w = getGeneration(ia, nabin, ahigh(0), atable(0, ie)) ! z direction, -1.0:upward, 0.0:horizontal, 1.0:downward
      sx = sqrt(1.0 - w**2) ! sin(theta)
      u = sx*cos(phi)  ! x direction
      v = sx*sin(phi)  ! y direction
      ener_out(i) = e
      u_out(i) = u
      v_out(i) = v
      w_out(i) = w
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end

function getGeneration(ibin, nbin, high, table)
   implicit real*8(a - h, o - z)
   dimension high(0:nbin)
   dimension table(0:nbin)

   rand = grnd() ! random number
   do i = 1, nbin - 1
      if (rand .le. table(i)) exit
   end do
   ibin = i ! bin ID

   rand = grnd() ! random number
   getGeneration = high(ibin - 1)*rand + high(ibin)*(1.0d0 - rand)

   return

end

! A C-program for MT19937: Real number version
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
subroutine sgrnd(seed)
!
   implicit integer(a - z)
!
! Period parameters
   parameter(N=624)
!
   dimension mt(0:N - 1)
!                     the array for the state vector
   common/block/mti, mt
   save/block/
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
   mt(0) = iand(seed, -1)
   do 1000 mti = 1, N - 1
      mt(mti) = iand(69069*mt(mti - 1), -1)
1000  continue
!
      return
   end
!***********************************************************************
   double precision function grnd()
!
      implicit integer(a - z)
!
! Period parameters
      parameter(N=624)
      parameter(N1=N + 1)
      parameter(M=397)
      parameter(MATA=-1727483681)
!                                    constant vector a
      parameter(UMASK=-2147483648)
!                                    most significant w-r bits
      parameter(LMASK=2147483647)
!                                    least significant r bits
! Tempering parameters
      parameter(TMASKB=-1658038656)
      parameter(TMASKC=-272236544)
!
      dimension mt(0:N - 1)
!                     the array for the state vector
      common/block/mti, mt
      save/block/
      data mti/N1/
!                     mti==N+1 means mt[N] is not initialized
!
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!                        mag01(x) = x * MATA for x=0,1
!
      TSHFTU(y) = ishft(y, -11)
      TSHFTS(y) = ishft(y, 7)
      TSHFTT(y) = ishft(y, 15)
      TSHFTL(y) = ishft(y, -18)
!
      if (mti .ge. N) then
!                       generate N words at one time
         if (mti .eq. N + 1) then
!                            if sgrnd() has not been called,
            call sgrnd(4357)
!                              a default initial seed is used
         end if
!
         do 1000 kk = 0, N - M - 1
            y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
            mt(kk) = ieor(ieor(mt(kk + M), ishft(y, -1)), mag01(iand(y, 1)))
1000        continue
            do 1100 kk = N - M, N - 2
               y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
               mt(kk) = ieor(ieor(mt(kk + (M - N)), ishft(y, -1)), mag01(iand(y, 1)))
1100           continue
               y = ior(iand(mt(N - 1), UMASK), iand(mt(0), LMASK))
               mt(N - 1) = ieor(ieor(mt(M - 1), ishft(y, -1)), mag01(iand(y, 1)))
               mti = 0
               end if
!
               y = mt(mti)
               mti = mti + 1
               y = ieor(y, TSHFTU(y))
               y = ieor(y, iand(TSHFTS(y), TMASKB))
               y = ieor(y, iand(TSHFTT(y), TMASKC))
               y = ieor(y, TSHFTL(y))
!
               if (y .lt. 0) then
                  grnd = (dble(y) + 2.0d0**32)/(2.0d0**32 - 1.0d0)
               else
                  grnd = dble(y)/(2.0d0**32 - 1.0d0)
               end if
!
               return
            end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine get_parma_particles_weights(type_list, ntypes, emin, emax, alt, glat, glong, weight_list)
               use :: iso_c_binding
               ! Calculate the weight of initial particles from PARMA at given altitude, energy range, date and position
               implicit none
               real(8) :: getr, getd, getspecangfinal, getspec, gethp
               integer npart, nebin, nabin
               parameter(npart=33) ! number of applicable particle
               parameter(nebin=1000) ! number of energy mesh (divided by log)
               parameter(nabin=100) ! number of angle mesh (divided by linear)
               real(8), intent(in) :: emin, emax
               real(8), intent(in) :: glat  ! Latitude (deg), -90 =< glat =< 90
               real(8), intent(in) :: glong ! Longitude (deg), -180 =< glong =< 180
               real(8), intent(in) :: alt
               integer, intent(in) :: ntypes ! number of different type of particles to include
               integer, intent(in), dimension(ntypes) :: type_list ! list of particles to includes
               real(8), intent(out), dimension(ntypes) :: weight_list
               real(8) :: alt_tmp
               real(8) :: ehigh(0:nebin), emid(nebin) ! higher and middle point of energy bin
               real(8) :: ahigh(0:nabin), amid(nabin) ! higher and middle point of angular bin
               real(8) :: etable(0:nebin)            ! probability table (0.0 for 0, 1.0 for nebin)
               real(8) :: atable(0:nabin, 0:nebin)    ! probability table (0.0 for 0, 1.0 for nabin)
               real(8) :: type_table(0:ntypes)            ! probability table
               !real(8) :: e_table(0:ntypes,0:naltbin,0:nebin)            ! probability table
               !real(8) :: a_table(0:ntypes,0:naltbin,0:nabin,0:nebin)    ! probability table
               real(8), allocatable :: a_table(:, :, :)
               real(8), allocatable :: e_table(:, :)
               real(8) :: summ_weights
integer :: IangPart(0:33) = (/1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 6/) ! Angular particle ID
               real(8) :: PPII = acos(-1.0)
               real(8) :: spec_val = 0.0, spec_angle = 0.0
               ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
               integer :: ip, nevent, i, i_alt, ia, ie, idummy, i_type
               real(8) :: g, amin, amax, amid_ia, emid_ie
               real(8) :: elog, estep, astep, alt_step, s, rr, d, e, w, sampled_alt
               real(8) :: getgeneration, getGeneration_index
               real(8) :: sum_g_e_p = 0., dummy_double = 0.
               integer :: index_type
               integer :: iyear     ! year
               integer :: imonth    ! Month
               integer :: iday      ! Date

               ALLOCATE (a_table(0:ntypes, 0:nabin, 0:nebin))
               ALLOCATE (e_table(0:ntypes, 0:nebin))

               !write (*, *), type_list
               !write (*, *), ntypes, emin, emax, alt, glat, glong

               iyear = 2020  ! Year
               imonth = 9    ! Month
               iday = 10      ! Date

               if (emin .lt. 0.01) then
                  write (*, *), "Energies below 1 MeV are not allowed."
                  call abort
               end if

               g = 0.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin

               ! Set Conditions (location, time, and local geometry)
               idummy = 0

               ! Set energy and angle ranges for generation
               amin = -1.0   ! Minimum cosine of particle
               amax = 1.0   ! Maximum cosine of particle

               ! Make energy and angle mesh and altitude mesh
               elog = log10(emin)
               estep = (log10(emax) - log10(emin))/nebin
               do ie = 0, nebin
                  ehigh(ie) = 10d0**elog
                  if (ie /= 0) emid(ie) = sqrt(ehigh(ie)*ehigh(ie - 1))
                  elog = elog + estep
               end do

               astep = (amax - amin)/nabin
               do ia = 0, nabin
                  ahigh(ia) = amin + astep*ia
                  if (ia /= 0) amid(ia) = (ahigh(ia) + ahigh(ia - 1))*0.5
               end do

               ! Make cumulative probability table for particle type

               a_table(:, :, :) = 0.0d0 ! initialization
               e_table(:, :) = 0.0d0
               type_table(:) = 0.0d0

               s = getHP(iyear, imonth, iday, idummy) ! Solar activity (W index) of the day
               rr = getr(glat, glong)               ! Vertical cut-off rigidity (GV)
               d = getd(alt, glat)  ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

               if (isnan(rr)) then
                  write (*, *), "ERROR : rr is nan. Aborting."
                  call abort
               end if
               if (isnan(s)) then
                  write (*, *), "ERROR : s is nan. Aborting."
                  call abort
               end if
               if (isnan(d)) then
                  write (*, *), "ERROR : d is nan. Aborting."
                  call abort
               end if

               do i_type = 1, ntypes

                  ip = type_list(i_type)

                  do ie = 1, nebin
                     do ia = 1, nabin
                        emid_ie = emid(ie)
                        spec_val = getSpec(ip, s, rr, d, emid_ie, g)
                        !write(*,*) spec_val
                        amid_ia = amid(ia)
                        if (amid_ia > 1.0) amid_ia = 0.999999999
                        if (amid_ia < -1.0) amid_ia = -0.999999999

                        spec_angle = getSpecAngFinal(iangpart(ip), s, rr, d, emid_ie, g, amid_ia)

                        if (isnan(spec_angle)) then
                           write (*, *) "spec_angle is nan. Aborting."
                           call abort
                        end if

                        if (isnan(spec_val)) then
                           write (*, *) "spec_val is nan. Aborting."
                           write (*, *), ip, iangpart(ip), s, rr, d, emid_ie, g, amid_ia, spec_val
                           call abort
                        end if

                  a_table(i_type, ia, ie) = a_table(i_type, ia - 1, ie) + spec_val*spec_angle*(2.0*PPII)*(ahigh(ia) - ahigh(ia - 1))

                        if (isnan(a_table(i_type, ia - 1, ie))) then
                           write (*, *) "one a_table value is nan, aborting."
                           call abort
                        end if

                     end do
                     e_table(i_type, ie) = e_table(i_type, ie - 1) + a_table(i_type, nabin, ie)*(ehigh(ie) - ehigh(ie - 1))

                  end do

                  type_table(i_type) = type_table(i_type - 1) + e_table(i_type, nebin)

                  weight_list(i_type) = e_table(i_type, nebin)

               end do

               ! normalisation to 1

               do i_type = 1, ntypes
                  type_table(i_type) = type_table(i_type)/type_table(ntypes)
               end do

               summ_weights = 0

               do i_type = 1, ntypes
                  summ_weights = summ_weights + weight_list(i_type)
               end do

               do i_type = 1, ntypes
                  weight_list(i_type) = weight_list(i_type)/summ_weights
               end do

            end subroutine get_parma_particles_weights

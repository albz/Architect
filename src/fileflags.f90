!*****************************************************************************************************!
!             Copyright 2014-2016 Alberto Marocchino, Francesco Massimo                               !
!*****************************************************************************************************!

!*****************************************************************************************************!
!  This file is part of architect.                                                                    !
!                                                                                                     !
!  Architect is free software: you can redistribute it and/or modify                                  !
!  it under the terms of the GNU General Public License as published by                               !
!  the Free Software Foundation, either version 3 of the License, or                                  !
!  (at your option) any later version.                                                                !
!                                                                                                     !
!  Architect is distributed in the hope that it will be useful,                                       !
!  but WITHOUT ANY WARRANTY; without even the implied warranty of                                     !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      !
!  GNU General Public License for more details.                                                       !
!                                                                                                     !
!  You should have received a copy of the GNU General Public License                                  !
!  along with architect.  If not, see <http://www.gnu.org/licenses/>.                                 !
!*****************************************************************************************************!

     subroutine getlun (nunit)
      logical lopen
      do i=350,399
         inquire(i,opened=lopen)
         if(.not.lopen) then
            nunit = i
            return
         endif
      enddo
      print *,'all units from 350 to 399 are open'
      stop
      end

!====================================================================
!   SET/UNSET FILE FLAGS
!====================================================================
	subroutine setFileFlag(fname)
	character(*),intent(in):: fname
	integer lun,getPID

	call getlun(lun)
	open(lun,file=trim(fname))
	if(trim(fname).eq.'__started__') then
		write(lun,*) getPID()
	endif
	close(lun,status='keep')

	end subroutine setFileFlag

!======================================================================
	subroutine unsetFileFlag(fname)
	character(*),intent(in):: fname
	integer lun

	call getlun(lun)
	open(lun,file=trim(fname))
	close(lun,status='delete')

	end subroutine unsetFileFlag

!======================================================================
	function getFileFlag(fname) result(exists)
	logical exists
	character(*) fname
!     inquire about file's existence:
	inquire (file = trim(fname), exist = exists)
	end function getFileFlag

!====================================================================
	function isNewRun() result(ans)
	logical ans,getFileFlag
	ans=.not.getFileFlag('==resume==')
	end function isNewRun
!======================================================================
	subroutine clearFileFlags

	end subroutine clearFileFlags

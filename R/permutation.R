# Copyright 2012 Paolo Martini <paolo.martini@unipd.it>
#
#
# This file is part of clipper.
#
# clipper is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License
# version 3 as published by the Free Software Foundation.
#
# clipper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with clipper. If not, see <http://www.gnu.org/licenses/>.

getperms <- function(n, nperms) {
  if (n > 50)
    total.perms <- Inf
  else
    total.perms=factorial(n)
  if(total.perms <= nperms){
    perms=permute(1:n)
    all.perms.flag=1
    nperms.act=total.perms
  }
  if(total.perms > nperms){
    perms=matrix(NA,nrow=nperms,ncol=n)
    for(i in 1:nperms){
      perms[i,]=sample(1:n, size=n)
    }
    all.perms.flag=0
    nperms.act=nperms
  }
  return(list(perms=perms, all.perms.flag=all.perms.flag, nperms.act=nperms.act))
}

permute<-function(elem) {
  if(length(elem) == 2)
    return(matrix(c(elem,elem[2],elem[1]),nrow=2))
  last.matrix<-permute(elem[-1])
  dim.last<-dim(last.matrix)
  new.matrix<-matrix(0,nrow=dim.last[1]*(dim.last[2]+1),ncol=dim.last[2]+1)
  for(row in 1:(dim.last[1])) {
    for(col in 1:(dim.last[2]+1))
      new.matrix[row+(col-1)*dim.last[1],]<-insertValue(last.matrix[row,],elem[1],col)
  }
  return(new.matrix)
}

insertValue<-function(vec,newval,pos) {
  if(pos == 1) return(c(newval,vec))
  lvec<-length(vec)
  if(pos > lvec) return(c(vec,newval))
  return(c(vec[1:pos-1],newval,vec[pos:lvec]))
}

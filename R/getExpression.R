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

setGeneric("getExpression", function(expr, classes) standardGeneric("getExpression"))

setMethod("getExpression",
          signature("matrix", "numeric"),
          function(expr, classes) {
            if (is.null(rownames(expr)))
              stop("Gene names not specified.")
            
            if (NCOL(expr) != length(classes))
              stop("Class vector length and sample number differs.")
            
            if (!all((classes == 1 | classes == 2)))
              stop("Class vector should be made by either '1' or '2'.")

            t(expr)
          })

setMethod("getExpression",
          signature("ExpressionSet", "numeric"),
          function(expr, classes) {
            getExpression(exprs(expr), classes)
          })


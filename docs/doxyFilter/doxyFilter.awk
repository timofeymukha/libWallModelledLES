#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
#    \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     doxyFilter.awk
#
# Description
#     Converts cocoon style sentinel strings into doxygen style strings
#
#     Assumes comment strings are formatted as follows
#         //- General description
#         //- with a continuation of general description
#         //  Detailed description
#         //  with more particulars
#     This should be re-formatted as the following
#         /*!
#          * \brief General description
#          * with a continuation of general description
#          *
#          * Detailed description
#          * with more particulars
#          */
#------------------------------------------------------------------------------

# state:  0=normal, 1=brief, 2=details
# indent: whitespace content prior to the opening "//-"
BEGIN {
    state = 0
    indent = ""
}

/^\s*\/\/-/ {
    if (state == 0)
    {
        # Changed from normal to brief (start of comment block)
        ## indent = substr($0, 1, index($0, "/")-1)
        indent = $0
        sub(/\S.*/, "", indent)
        printf indent "/*!\n"
        printf indent " * \\brief "
        sub(/^\s*\/\/-\s*/, "")
        state = 1
    }
    else
    {
        # Within brief: replace leading space with proper indent amount
        printf indent
        sub(/^\s*\/\/-\s*/, " * ")
    }

    print
    next
}


/^\s*\/\// {
    if (state == 1)
    {
        # Change from brief to details. Extra line to start new paragraph.
        printf indent " *\n"
        state = 2
    }

    if (state == 2)
    {
        # Within details
        printf indent

        # '//' with 4 spaces or more - assume indent is intentional
        if (match($0, /^\s*\/\/(    )+/))
        {
            sub(/^\s*\/\/\s/, " *")
        }
        else
        {
            sub(/^\s*\/\/\s*/, " * ")
        }
    }

    print
    next
}


{
    # End comment filtering
    if (state)
    {
        printf indent " */\n"
        state = 0
    }
    print
    next
}

#------------------------------------------------------------------------------

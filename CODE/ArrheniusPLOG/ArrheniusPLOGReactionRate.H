/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ArrheniusPLOGReactionRate

Description
    PLOG in OpenFOAM.

    Usage:
        ArrheniusPLOGReactionRateCoeffs
        {
            A         9.2E+38; // place holder
            beta      -9.43;   // place holder
            Ta        7000.;   // place holder
            ArrheniusData      // optional
            (
                (0.1    9.2E+38 -8.65 7000.  )   // PLOG /p (Pa) A beta Ta/
                (1.0    9.5E+42 -9.43 11200. )   // PLOG /p A beta Ta/
                (10.    1.5E+45 -9.69 15100. )   // PLOG /p A beta Ta/
                (100.   1.8E+43 -8.78 16800. )   // PLOG /p A beta Ta/
                (1.0E+5 4.4E+09  1.45 2400.  )   // PLOG /p A beta Ta/
            );
        }
        
SourceFiles
    ArrheniusPLOGReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusPLOGReactionRate_H
#define ArrheniusPLOGReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"
#include "FixedList.H"

#include "SortableList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class ArrheniusPLOGReactionRate;

Ostream& operator<<(Ostream&, const ArrheniusPLOGReactionRate&);


/*---------------------------------------------------------------------------*\
                   Class ArrheniusPLOGReactionRate Declaration
\*---------------------------------------------------------------------------*/

class ArrheniusPLOGReactionRate
{
    // Private Data

        // place holder
        scalar A_;
        scalar beta_;
        scalar Ta_;
        

    // Helper class to store specie-local volatile data
    class ArrheniusData
    {
        // Private Data

            //- pressure range
            scalar p_;

            scalar A_;
            scalar beta_;
            scalar Ta_;

    public:

        // Constructors

            //- Null constructor
            ArrheniusData()
            :
                p_(0.0),
                A_(0.0),
                beta_(0.0),
                Ta_(0.0)
            {}

            //- Construct from Istream
            ArrheniusData(Istream& is)
            :
                p_(readScalar(is)),
                A_(readScalar(is)),
                beta_(readScalar(is)),
                Ta_(readScalar(is))
            {}
/*
            //- Construct from a copy of ArrheniusData
            ArrheniusData(ArrheniusData& pcm)
            :
                p_(pcm.p()),
                A_(pcm.A()),
                beta_(pcm.beta()),
                Ta_(pcm.Ta())
            {}
*/
        //- Destructor
        ~ArrheniusData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const scalar& p() const
                {
                    return p_;
                }

                //- Return const access to the pre-exponential factor
                scalar A() const
                {
                    return A_;
                }

                //- Return const access to beta
                scalar beta() const
                {
                    return beta_;
                }

                //- Return const access to the activation energy
                scalar Ta() const
                {
                    return Ta_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, ArrheniusData& ad)
            {
                is.readBeginList("ArrheniusData");
                is  >>  ad.p_  >> ad.A_  >> ad.beta_ >> ad.Ta_;
                is.readEndList("ArrheniusData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const ArrheniusData& ad)
            {
                os  << token::BEGIN_LIST
                    << ad.p_ << token::SPACE
                    << ad.A_ << token::SPACE
                    << ad.beta_ << token::SPACE
                    << ad.Ta_
                    << token::END_LIST;

                return os;
            }
    };

    // Private Data

        // Model constants
            //- List of ArrheniusData data - (p A bata Ta)
            List<ArrheniusData> ArrheniusData_;

public:

    // Constructors

        //- Construct from components
        inline ArrheniusPLOGReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const List<ArrheniusData> ads
        );

        //- Construct from dictionary
        inline ArrheniusPLOGReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "ArrheniusPLOG";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        inline const List<Tuple2<label, scalar>>& beta() const;

        //- Species concentration derivative of the pressure dependent term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const ArrheniusPLOGReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ArrheniusPLOGReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

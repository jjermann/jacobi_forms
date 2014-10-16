r"""
Modular forms for Hecke triangle groups

AUTHORS:

- Jonas Jermann (2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ, infinity

from sage.modules.module import Module
from sage.categories.all import Modules
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from hecke_triangle_groups import HeckeTriangleGroup
from abstract_space import FormsSpace_abstract

def canonical_parameters(group, base_ring, k, ep, n=None):
    r"""
    Return a canonical version of the parameters.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.space import canonical_parameters
        sage: canonical_parameters(5, ZZ, 20/3, int(1))
        (Hecke triangle group for n = 5, Integer Ring, 20/3, 1, 5)

        sage: canonical_parameters(infinity, ZZ, 2, int(-1))
        (Hecke triangle group for n = +Infinity, Integer Ring, 2, -1, +Infinity)
    """

    if not (n is None):
        group = n

    if (group == infinity):
        group = HeckeTriangleGroup(infinity)
    else:
        try:
            group = HeckeTriangleGroup(ZZ(group))
        except TypeError:
            group = HeckeTriangleGroup(group.n())

    n = group.n()
    k = QQ(k)
    if (ep == None):
        if (n == infinity):
            ep = (-1)**(k/ZZ(2))
        elif (ZZ(2).divides(n)):
            ep = (-1)**(k*ZZ(n-2)/ZZ(4))
        else:
            ep = (-1)**(k*ZZ(n-2)/ZZ(2))
    ep = ZZ(ep)

    if (n == infinity):
        num = (k-(1-ep)) / ZZ(4)
    else:
        num = (k-(1-ep)*n/(n-2)) * (n-2) / ZZ(4)

    try:
        num = ZZ(num)
    except TypeError:
        pass
        #raise ValueError("Invalid or non-occuring weight k={}, ep={}!".format(k,ep))

    return (group, base_ring, k, ep, n)


class QuasiMeromorphicJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi meromorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiMeromorphicModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, ZZ, 20/3, int(1))
            sage: QuasiMeromorphicModularForms(5, ZZ, 20/3, int(1)) == QuasiMeromorphicModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi meromorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MF = QuasiMeromorphicModularForms(5, ZZ, 20/3, 1)
            sage: MF
            QuasiMeromorphicModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi meromorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_space() == MF
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "quasi", "mero"])

class QuasiWeakJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi weakly holomorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiWeakModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(4, ZZ, 8, -1)
            sage: QuasiWeakModularForms(4, ZZ, 8, -1) == QuasiWeakModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi weakly holomorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiWeakModularForms
            sage: MF = QuasiWeakModularForms(4, ZZ, 8, 1)
            sage: MF
            QuasiWeakModularForms(n=4, k=8, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi weakly holomorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "quasi", "weak"])

class QuasiJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, ZZ, 10/3, -1)
            sage: QuasiModularForms(5, ZZ, 10/3) == QuasiModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(5, ZZ, 20/3, 1)
            sage: MF
            QuasiModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "quasi", "holo"])

class QuasiCuspJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi cusp forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiCuspForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(8, ZZ, 16/3, None)
            sage: QuasiCuspForms(8, ZZ, 16/3) == QuasiCuspForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi cusp forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(8, ZZ, 16/3)
            sage: MF
            QuasiCuspForms(n=8, k=16/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi cuspidal
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True

            sage: QuasiCuspForms(n=infinity)
            QuasiCuspForms(n=+Infinity, k=0, ep=1) over Integer Ring
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "quasi", "cusp"])

class MeromorphicJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) meromorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, MeromorphicModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(3, ZZ, 0, 1)
            sage: MeromorphicModularForms() == MeromorphicModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) meromorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import MeromorphicModularForms
            sage: MF = MeromorphicModularForms()
            sage: MF
            MeromorphicModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF.analytic_type()
            meromorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "mero"])

class WeakJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) weakly holomorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, WeakModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, CC, 20/3, None)
            sage: WeakModularForms(5, CC, 20/3) == WeakModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) weakly holomorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(5, CC, 20/3)
            sage: MF
            WeakModularForms(n=5, k=20/3, ep=1) over Complex Field with 53 bits of precision
            sage: MF.analytic_type()
            weakly holomorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "weak"])

class JacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, ModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(3, ZZ, 0, None)
            sage: ModularForms() == ModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms()
            sage: MF
            ModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF.analytic_type()
            modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.module()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True

            sage: MF = ModularForms(n=infinity, k=8)
            sage: MF
            ModularForms(n=+Infinity, k=8, ep=1) over Integer Ring
            sage: MF.analytic_type()
            modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type = self.AT(["jacobi", "holo"])

class CuspJacobiForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) cusp forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, CuspForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(6, ZZ, 6, 1)
            sage: CuspForms(6, ZZ, 6, 1) == CuspForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) cusp forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(6, ZZ, 6, 1)
            sage: MF
            CuspForms(n=6, k=6, ep=1) over Integer Ring
            sage: MF.analytic_type()
            cuspidal
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.module()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["jacobi", "cusp"])

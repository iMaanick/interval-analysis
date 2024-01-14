import math as m
from numbers import Number


class Interval:
    def __init__(self, left: int | float, right: int | float) -> 'Interval':
        self.__left = left
        self.__right = right

    def __str__(self) -> str:
        return '[' + str(self.__left) + ', ' + str(self.__right) + ']'

    def __repr__(self) -> str:
        return '[' + str(self.__left) + ', ' + str(self.__right) + ']'

    def __add__(self, other: 'Interval') -> 'Interval':
        return Interval(self.__left + other.__left, self.__right + other.__right)

    def __radd__(self, other: 'Interval') -> 'Interval':
        return Interval(self.__left + other.__left, self.__right + other.__right)

    def __sub__(self, other: 'Interval') -> 'Interval':
        return Interval(self.__left - other.__right, self.__right - other.__left)

    def __rsub__(self, other: 'Interval') -> 'Interval':
        return Interval(self.__left - other.__right, self.__right - other.__left)

    def __mul__(self, other: 'Interval') -> 'Interval':
        return Interval(min(self.__left * other.__left, self.__left * other.__right,
                            self.__right * other.__left, self.__right * other.__right),
                        max(self.__left * other.__left, self.__left * other.__right,
                            self.__right * other.__left, self.__right * other.__right))

    def __rmul__(self, other: 'Interval') -> 'Interval':
        return Interval(min(self.__left * other.__left, self.__left * other.__right,
                            self.__right * other.__left, self.__right * other.__right),
                        max(self.__left * other.__left, self.__left * other.__right,
                            self.__right * other.__left, self.__right * other.__right))

    def __truediv__(self, other: 'Interval') -> 'Interval':
        if any([m.isclose(other.__left, 0.0), m.isclose(other.__left, 0.0),
                m.isclose(other.__right, 0.0), m.isclose(other.__right, 0.0)]):
            raise ZeroDivisionError("Divide by zero.")

        return Interval(min(self.__left / other.__left, self.__left / other.__right,
                            self.__right / other.__left, self.__right / other.__right),
                        max(self.__left / other.__left, self.__left / other.__right,
                            self.__right / other.__left, self.__right / other.__right))

    def __rtruediv__(self, other: 'Interval') -> 'Interval':
        if any([m.isclose(other.__left, 0.0), m.isclose(other.__left, 0.0),
                m.isclose(other.__right, 0.0), m.isclose(other.__right, 0.0)]):
            raise ZeroDivisionError("Divide by zero.")

        return Interval(min(self.__left / other.__left, self.__left / other.__right,
                            self.__right / other.__left, self.__right / other.__right),
                        max(self.__left / other.__left, self.__left / other.__right,
                            self.__right / other.__left, self.__right / other.__right))

    def __contains__(self, num: int | float) -> bool:
        return self.__left <= num and num <= self.__right

    @property
    def left(self) -> int | float:
        return self.__left

    @left.setter
    def left(self, num: int | float) -> None:
        if not isinstance(num, Number):
            raise ValueError("num must be numeric.")
        self.__left = num

    @left.deleter
    def left(self) -> None:
        raise PermissionError("left is not deletable")

    @property
    def right(self) -> int | float:
        return self.__right

    @right.setter
    def right(self, num: int | float) -> None:
        if not isinstance(num, Number):
            raise ValueError("num must be numeric.")
        self.__right = num

    @right.deleter
    def right(self) -> None:
        raise PermissionError("right is not deletable")

    @property
    def mid(self) -> int | float:
        return (self.__left + self.__right) / 2

    @mid.setter
    def mid(self, *args) -> None:
        raise PermissionError("mid is not editable")

    @mid.deleter
    def mid(self) -> None:
        raise PermissionError("mid is not deletable")

    @property
    def rad(self) -> int | float:
        return (self.__right - self.__left) / 2

    @rad.setter
    def rad(self, *args) -> None:
        raise PermissionError("rad is not editable")

    @rad.deleter
    def rad(self) -> None:
        raise PermissionError("rad is not deletable")

    @property
    def wid(self) -> int | float:
        return self.__right - self.__left

    @wid.setter
    def wid(self, *args) -> None:
        raise PermissionError("wid is not editable")

    @wid.deleter
    def wid(self) -> None:
        raise PermissionError("wid is not deletable")


empty = Interval(0, 0)
infinite = Interval(-m.inf, m.inf)
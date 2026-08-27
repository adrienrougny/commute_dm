"""Register pylpg node classes for the `momapy_bel` element classes.

Imported for its effect alone, wherever a `BELModel` is read back out of the
database: `execute_query_as_objects` can only rebuild an object whose momapy
class has a registered node class, and a session registers only what it saved
or what a static type hint names. This is a copy of
`momapy_kb.lpg.celldesigner` pointed at `momapy_bel.core`, BEL having no such
module of its own (and no layout).
"""

import sys

import fieldz_kb.lpg.core

import momapy.core.map
import momapy.core.model
import momapy.core.elements

import momapy_bel.core

import momapy_kb.lpg.types

ctx = fieldz_kb.lpg.core.get_default_context()
momapy_kb.lpg.types.register_momapy_plugins(ctx)
module = momapy_bel.core
for attr_name in dir(module):
    if not attr_name.startswith("_"):
        attr_value = getattr(module, attr_name)
        if isinstance(attr_value, type) and issubclass(
            attr_value,
            (
                momapy.core.elements.ModelElement,
                momapy.core.map.Map,
                momapy.core.model.Model,
            ),
        ):
            node_class = fieldz_kb.lpg.core.get_or_make_node_class_from_type(
                ctx, attr_value, make_node_classes_recursively=True
            )
            setattr(sys.modules[__name__], node_class.__name__, node_class)

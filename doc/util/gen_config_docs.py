'''
Generate documentation from a pydantic model's type hints and field docstrings
in yaml or json format. Requires 'pydantic'.

Example:
    # foo.py
    from pydantic import BaseModel

    class Bar(BaseModel):
        baz: int
        """
        The baz field holds the number of quox in the observable universe.
        """

    # from command line
    python gen_config_docs.py foo.Bar
    # The baz field holds the number of quox in the observable universe.
    baz: int
'''

from __future__ import annotations

import argparse
import ast
import dataclasses
import inspect
import itertools
import typing

import pydantic


@dataclasses.dataclass
class Field:
    name: str
    annotation: str
    default: str | None = None
    docstring: str | None = None


@dataclasses.dataclass
class PydanticField:
    field: Field
    fields: list[Field | PydanticField] = dataclasses.field(default_factory=list)


def assert_never(arg: typing.Never):
    raise AssertionError(f"Expected code to be unreachable, but got: {arg}")


def format_attribute(attr: ast.Attribute) -> str:
    if type(attr.value) == ast.Name:
        return f"{attr.value.id}.{attr.attr}"
    return attr.attr


def format_fn(fn: ast.Call) -> str:
    name = format_annotation(fn.func)
    args = [format_annotation(arg) for arg in fn.args]
    kwargs = [f"{kw.arg}={format_annotation(kw.value)}" for kw in fn.keywords]
    return f"{name}({', '.join(itertools.chain(args, kwargs))})"


def format_unary(val: ast.UnaryOp) -> str:
    if type(val.op) == ast.UAdd:
        op = "+"
    elif type(val.op) == ast.USub:
        op = "-"
    else:
        assert_never(type(val))
    assert type(val.operand) == ast.Constant
    return f"{op}{val.operand.value}"


def format_annotation(
    sub: (
        ast.Name
        | ast.Attribute
        | ast.Constant
        | ast.Tuple
        | ast.Call
        | ast.UnaryOp
        | ast.Lambda
    ),
) -> str:
    if type(sub) == ast.Name:
        return sub.id
    elif type(sub) == ast.Attribute:
        return format_attribute(sub)
    elif type(sub) == ast.Constant:
        return repr(sub.value)
    elif type(sub) == ast.UnaryOp:
        return format_unary(sub)
    elif type(sub) == ast.Subscript:
        return format_subscript(sub)
    elif type(sub) == ast.Call:
        return format_fn(sub)
    if type(sub) == ast.Lambda:
        # TODO: this could be improved
        return "lambda"
    elif type(sub) == ast.Tuple:
        return f"{', '.join([format_annotation(el) for el in sub.elts])}"

    import astpretty

    astpretty.pprint(sub)
    assert_never(type(sub))


def format_subscript(sub: ast.Subscript) -> str:
    if type(sub.value) == ast.Name:
        attr = sub.value.id
    elif type(sub.value) == ast.Attribute:
        attr = format_attribute(sub.value)
    elif type(sub.value) == ast.Constant:
        attr = sub.value
    else:
        assert_never(sub.value)
    return f"{attr}[{format_annotation(sub.slice)}]"


def is_pydantic_subtype(node: ast.ClassDef) -> bool:
    def predicate(node: ast.Name | ast.Attribute | ast.expr) -> bool:
        if type(node) == ast.Name and "BaseModel" in node.id:
            return True
        elif type(node) == ast.Attribute and "BaseModel" in node.attr:
            return True
        return False

    return any([b for b in node.bases if predicate(b)])


class NodeVisitor(ast.NodeVisitor):
    def __init__(self) -> None:
        self.models: dict[str, list[Field]] = {}

    def visit_ClassDef(self, node: ast.ClassDef) -> typing.Any:
        if is_pydantic_subtype(node):
            fields = []
            for i, item in enumerate(node.body):
                item = node.body[i]
                if type(item) != ast.AnnAssign:
                    continue
                name = item.target.id
                ann = format_annotation(item.annotation)
                default = format_annotation(item.value) if item.value else None
                if i + 1 < len(node.body) and type(node.body[i + 1]) == ast.Expr:
                    doc = node.body[i + 1].value.value
                else:
                    doc = None
                fields.append(Field(name, ann, default, doc))
            self.models[node.name] = fields

        self.generic_visit(node)


T = typing.TypeVar("T", covariant=True)


def first(l: list[T], pred: typing.Callable[[T], bool]) -> int:
    for i, item in enumerate(l):
        if pred(item):
            return i
    return -1


def get_model_fields(
    model: pydantic.BaseModel,
) -> list[Field | PydanticField]:
    nv = NodeVisitor()
    try:
        nv.visit(ast.parse(inspect.getsource(model)))
    except AssertionError as e:
        print(f"model: {model.__module__}.{model.__name__}")
        raise e
    fields = nv.models[model.__name__]
    for field in model.__fields__.values():
        # TODO: test with future annotations
        type_ = field.type_ if hasattr(field, "type_") else field.annotation
        if inspect.isclass(type_) and issubclass(type_, pydantic.BaseModel):
            idx = first(fields, lambda f: isinstance(f, Field) and f.name == field.name)
            assert idx != -1
            subfields = get_model_fields(field.type_)
            fields[idx] = PydanticField(fields[idx], subfields)
    return fields


def pretty_print(l: list[Field | PydanticField], indent: str = ""):
    for item in l:
        if isinstance(item, Field):
            print(f"{indent}{item}")
        else:
            print(f"{indent}{item.field.name}: {item.field.annotation}")
            pretty_print(item.fields, indent + " " * 2)


def json_print(l: list[Field | PydanticField], indent: str = ""):
    for item in l:
        if isinstance(item, Field):
            if item.docstring is not None:
                docstring = item.docstring.strip().replace("\n", "\n# ")
                print(f"{indent}# {docstring}")
            print(f"{indent}{item.name}: {item.annotation}")
        else:
            if item.field.docstring is not None:
                docstring = item.field.docstring.strip().replace("\n", "\n# ")
                print(f"{indent}# {docstring}")
            print(f"{indent}{item.field.name}: {item.field.annotation}")
            json_print(item.fields, indent + " " * 2)


class DictField(typing.TypedDict):
    type: str
    docs: str | None


class DictFields(typing.TypedDict):
    type: str
    docs: str | None
    fields: dict[str, DictField | DictFields]


def field_as_dict(field: Field) -> dict[str, DictField]:
    return {
        field.name: {
            "type": field.annotation,
            "docs": field.docstring,
        }
    }


def pydantic_field_as_dict(field: PydanticField) -> dict[str, DictFields]:
    fields = {}
    for f in field.fields:
        fields.update(as_dict(f))

    return {
        field.field.name: {
            "type": field.field.annotation,
            "docs": field.field.docstring,
            "fields": fields,
        }
    }


def as_dict(
    field: Field | PydanticField,
) -> dict[str, DictField] | dict[str, DictFields]:
    if type(field) == Field:
        return field_as_dict(field)
    elif type(field) == PydanticField:
        return pydantic_field_as_dict(field)

    assert_never(type(field))


def print_json(fields: list[PydanticField | Field]):
    dict_fields = {}
    for f in fields:
        dict_fields.update(as_dict(f))

    import json

    print(json.dumps(dict_fields, indent=2))


def print_ast(o: object) -> None:
    tree = ast.parse(inspect.getsource(o))
    import astpretty

    astpretty.pprint(tree)


def _print_yaml_docstring(docstring: str, indent: str):
    docstring = docstring.strip().replace("\n", "\n# ")
    print(f"{indent}# {docstring}")


def _print_yaml_field(field: Field, indent: str):
    default = f" = {field.default}" if field.default else ""
    print(f"{indent}{field.name}: {field.annotation}{default}")


def print_yaml(l: list[Field | PydanticField], indent: str = ""):
    for item in l:
        if isinstance(item, Field):
            if item.docstring is not None:
                _print_yaml_docstring(item.docstring, indent)
            _print_yaml_field(item, indent)
        else:
            if item.field.docstring is not None:
                _print_yaml_docstring(item.field.docstring, indent)
            _print_yaml_field(item.field, indent)
            print_yaml(item.fields, indent + " " * 2)


def import_model(model_name: str) -> pydantic.BaseModel:
    import importlib

    module, model = model_name.rsplit(".", 1)
    mod = importlib.import_module(module)
    model: pydantic.BaseModel = getattr(mod, model)

    assert issubclass(model, pydantic.BaseModel)
    return model


def main() -> int:
    description = "Generate a pydantic model's documentation using its type hints and field docstrings"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "model", nargs=1, help="python formatted import path to pydantic model"
    )
    parser.add_argument(
        "--ast",
        action="store_true",
        help="output abstract syntax tree (requires 'astpretty' library)",
    )
    parser.add_argument("--json", action="store_true", help="use json format")
    args = parser.parse_args()

    model = import_model(args.model[0])
    if args.ast:
        print_ast(model)
        return 0

    fields = get_model_fields(model)
    if args.json:
        print_json(fields)
    else:
        print_yaml(fields)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

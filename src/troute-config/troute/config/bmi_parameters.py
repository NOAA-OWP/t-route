from pydantic import BaseModel


class BMIParameters(BaseModel):
    segment_number: int
    waterbody_number: int
    io_number: int
    upstream_number: int

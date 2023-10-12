from pydantic import BaseModel, Field
from typing import Optional, List


class BMIParameters(BaseModel, extra='forbid'):
    flowpath_columns: Optional[List[str]] = Field(
        default_factory=lambda: [
            'id',
            'toid',
            'lengthkm',
            ]
        )
    attributes_columns: Optional[List[str]] = Field(
        default_factory=lambda: [
            'attributes_id',
            'rl_gages',
            'rl_NHDWaterbodyComID',
            'MusK',
            'MusX',
            'n',
            'So',
            'ChSlp',
            'BtmWdth',
            'nCC',
            'TopWdthCC',
            'TopWdth',
            ]
        )
    waterbody_columns: Optional[List[str]] = Field(
        default_factory=lambda: [
            'hl_link',
            'ifd',
            'LkArea',
            'LkMxE',
            'OrificeA',
            'OrificeC',
            'OrificeE',
            'WeirC',
            'WeirE',
            'WeirL',
            ]
        )
    network_columns: Optional[List[str]] = Field(
        default_factory=lambda: [
            'network_id',
            'hydroseq',
            'hl_uri',
            ]
        )
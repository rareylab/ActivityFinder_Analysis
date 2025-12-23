from enum import Enum


class SequenceTypes(Enum):
    seqres = 1
    alignedseqres = 2
    alignedseqresafterbm = 3
    seqresseqadv = 4
    alignedseqresseqadv = 5
    alignedseqresseqadvafterbm = 6
    atomseq = 7
    alignedatomseq = 8
    alignedatomseqafterbm = 9
    alignedatomseqseqadv = 10
    alignedatomseqseqadvafterbm = 11
    chembl = 12
    uniprot = 13
    trembl = 14
    manual = 15
    varseq = 16
    compseq = 17

    def __new__(cls, value):
        member = object.__new__(cls)
        member._value_ = value
        return member

    def __int__(self):
        return self.value


class StringLookUp:
    type = {
        SequenceTypes.seqres: "seqres",
        SequenceTypes.seqresseqadv: "seqresseqadv",
        SequenceTypes.atomseq: "atomseq",
        SequenceTypes.chembl: "chembl",
        SequenceTypes.uniprot: "uniprot",
        SequenceTypes.trembl: "trembl",
        SequenceTypes.manual: "manual",
    }

    string_sequence_type = {
        "seqres": SequenceTypes.seqres,
        "seqresseqadv": SequenceTypes.seqresseqadv,
        "atomseq": SequenceTypes.atomseq,
        "chembl": SequenceTypes.chembl,
        "uniprot": SequenceTypes.uniprot,
        "trembl": SequenceTypes.trembl,
        "manual": SequenceTypes.manual,
    }

    @classmethod
    def get_string_for_sequence_type(cls, sequence_type: str):
        if sequence_type in cls.type.values():
            return cls.string_sequence_type[sequence_type]
        raise KeyError(
            f"sequence type {sequence_type} unknown. "
            f'available: ({",".join([type for type in cls.type.values()])})'
        )

import xml.etree.ElementTree as ET

xml = """
<proteome xmlns="http://uniprot.org/proteome" proteinCount="2778">
<upid>UP000722097</upid>
<taxonomy>2030927</taxonomy>
<modified>2023-05-02</modified>
<isReferenceProteome>false</isReferenceProteome>
<isRepresentativeProteome>false</isRepresentativeProteome>
<genomeAssembly>
<genomeAssemblySource>ENA/EMBL</genomeAssemblySource>
<genomeAssembly>GCA_012729335.1</genomeAssembly>
<genomeAssemblyUrl>https://www.ebi.ac.uk/ena/browser/view/GCA_012729335.1</genomeAssemblyUrl>
<genomeRepresentation>full</genomeRepresentation>
</genomeAssembly>
<genomeAnnotation>
<genomeAnnotationSource>ENA/EMBL</genomeAnnotationSource>
<genomeAnnotationUrl>https://www.ebi.ac.uk/ena/browser/view/GCA_012729335.1</genomeAnnotationUrl>
</genomeAnnotation>
<component name="Unassembled WGS sequence" proteinCount="2778">
<description>Bacteroidales bacterium</description>
<biosampleId>SAMN13893707</biosampleId>
<genomeAccession>JAAYDJ010000000</genomeAccession>
<genomeAnnotation>
<genomeAnnotationSource>ENA/EMBL</genomeAnnotationSource>
</genomeAnnotation>
</component>
<annotationScore normalizedAnnotationScore="2"/>
<scores name="busco">
<property name="completed" value="509"/>
<property name="completedSingle" value="505"/>
<property name="completedDuplicated" value="4"/>
<property name="total" value="541"/>
<property name="fragmented" value="3"/>
<property name="missing" value="29"/>
<property name="score" value="94"/>
<property name="lineage" value="bacteroidales_odb10"/>
</scores>
<scores name="cpd">
<property name="averageCds" value="0"/>
<property name="confidence" value="0"/>
<property name="proteomeCount" value="0"/>
<property name="stdCdss" value="0.0"/>
<property name="status" value="Unknown"/>
</scores>
</proteome>
"""

root = ET.fromstring(xml)
details = {
    'proteinCount': root.attrib['proteinCount'],
    'upid': root.find('{http://uniprot.org/proteome}upid').text,
    'taxonomy': root.find('{http://uniprot.org/proteome}taxonomy').text,
    'modified': root.find('{http://uniprot.org/proteome}modified').text,
    'isReferenceProteome': root.find('{http://uniprot.org/proteome}isReferenceProteome').text,
    'isRepresentativeProteome': root.find('{http://uniprot.org/proteome}isRepresentativeProteome').text
}
for scores in root.findall('{http://uniprot.org/proteome}scores'):
    score_name = scores.attrib['name']
    for property in scores.findall('{http://uniprot.org/proteome}property'):
        prop_name = property.attrib['name']
        prop_value = property.attrib['value']
        details[f'{score_name}_{prop_name}'] = prop_value

print(details)
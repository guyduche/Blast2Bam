<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>
<xsl:param name="fileType"/>

<xsl:variable name="mit2">/*
The MIT License (MIT)

Copyright (c) 2015 Aurelien Guy-Duche

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2015 creation

*/</xsl:variable>

<xsl:template match="/">
<xsl:choose>
<xsl:when test="$fileType='h'">
<xsl:value-of select="$mit2"/>
#ifndef PARSEXML_H_INCLUDED
#define PARSEXML_H_INCLUDED

#include &lt;libxml/xmlreader.h&gt;
</xsl:when>
<xsl:when test="$fileType='c'">
<xsl:value-of select="$mit2"/>
#include &lt;string.h&gt;
#include &lt;stdlib.h&gt;
#include "parseXML.h"
#include "utils.h"

#define PARSEMACRO(TAG, STRUCT, FIELD, FUNCT) (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*)TAG)) \
        { \
            if (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*)"blastOutput_iterations")) \
                end = 1; /* Useful to stop parseBlastOutput() before the first iteration of results */ \
            else \
            { \
                evt = safeXmlTextReaderRead(reader); \
                STRUCT->FIELD = FUNCT; \
                evt = safeXmlTextReaderRead(reader); \
            } \
        }

</xsl:when>
</xsl:choose>

<xsl:apply-templates select="schema"/>

</xsl:template>

<xsl:template match="schema">
<xsl:apply-templates select="struct" mode="definition"/>
<xsl:if test="$fileType='h'">
/************************************************************************************/
/*  Prototypes                                                                      */
/************************************************************************************/
<xsl:apply-templates select="struct" mode="prototype"/>
#endif // PARSEXML_H_INCLUDED
</xsl:if>
</xsl:template>

<xsl:template match="struct" mode="definition">
<xsl:if test="$fileType='h'">
/* start define <xsl:value-of select="@name"/> */
typedef struct <xsl:value-of select="@name"/>
{
    <xsl:for-each select="field">
    <xsl:choose>
        <xsl:when test="@type='string'">char*</xsl:when>
        <xsl:when test="@type='int'">int</xsl:when>
        <xsl:when test="@type='double'">double</xsl:when>
        <xsl:otherwise><xsl:value-of select="@type"/>*</xsl:otherwise>		
    </xsl:choose><xsl:text> </xsl:text><xsl:apply-templates select="." mode="cname"/>;
    </xsl:for-each>
    struct <xsl:value-of select="@name"/>* next;
} <xsl:value-of select="@name"/>, *<xsl:value-of select="@name"/>Ptr;
</xsl:if>

<xsl:if test="$fileType='c'">
<xsl:if test="@name!='BlastOutput' and @name!='Iteration'">static<xsl:text> </xsl:text></xsl:if><xsl:value-of select="@name"/>Ptr parse<xsl:value-of select="@name"/>(xmlTextReaderPtr reader)
{
    int evt = 1, end = 0;
    <xsl:value-of select="@name"/>Ptr p = (<xsl:value-of select="@name"/>Ptr) safeCalloc(1, sizeof(<xsl:value-of select="@name"/>));

    evt = safeXmlTextReaderRead(reader);

    if (xmlTextReaderNodeType(reader) == XML_READER_TYPE_END_ELEMENT)
        end = 1;
    else
        evt = safeXmlTextReaderRead(reader);

    while (evt == 1 &amp;&amp; !end)
    {	
        if (xmlTextReaderNodeType(reader) == XML_READER_TYPE_END_ELEMENT)
            end = 1;
        <xsl:for-each select="field">
        else if PARSEMACRO("<xsl:value-of select="@name"/>", p, <xsl:apply-templates select="." mode="cname"/>, <xsl:choose>
                <xsl:when test="@type='string'">safeStrdup((char*)xmlTextReaderConstValue(reader))</xsl:when>
                <xsl:when test="@type='int'">(int) strtol((char*)xmlTextReaderConstValue(reader), NULL, 10)</xsl:when>
                <xsl:when test="@type='double'">strtod((char*)xmlTextReaderConstValue(reader), NULL)</xsl:when>
                <xsl:otherwise>parse<xsl:value-of select="@type"/>(reader)</xsl:otherwise>
            </xsl:choose>)
        </xsl:for-each>

        evt = safeXmlTextReaderRead(reader);
    }

    if (xmlTextReaderNodeType(reader) != XML_READER_TYPE_END_ELEMENT)
        evt = safeXmlTextReaderRead(reader);
    <xsl:if test="@name!='Iteration'">
    if (!xmlStrcasecmp(xmlTextReaderConstName(reader), (xmlChar*) "<xsl:value-of select="@name"/>"))
        p->next = parse<xsl:value-of select="@name"/>(reader);
    </xsl:if>
    return p;
}

<xsl:if test="@name!='BlastOutput' and @name!='Iteration'">static<xsl:text> </xsl:text></xsl:if>void dealloc<xsl:value-of select="@name"/>(<xsl:value-of select="@name"/>Ptr element)
{
    <xsl:if test="@name='BlastOutput' or @name='Iteration' or @name='Hit'">void* prec = NULL;</xsl:if>
    <xsl:for-each select="field">
    <xsl:choose>
        <xsl:when test="@type='int'"></xsl:when>
        <xsl:when test="@type='double'"></xsl:when>
        <xsl:when test="@type='string'">
    if (element-><xsl:apply-templates select="." mode="cname"/> != NULL) free(element-><xsl:apply-templates select="." mode="cname"/>);
        </xsl:when>
        <xsl:otherwise>
    while(element-><xsl:apply-templates select="." mode="cname"/> != NULL)
    {
        prec = element-><xsl:apply-templates select="." mode="cname"/>;
        element-><xsl:apply-templates select="." mode="cname"/> = element-><xsl:apply-templates select="." mode="cname"/>->next;
        dealloc<xsl:value-of select="@type"/>(prec);
    }
        </xsl:otherwise>
    </xsl:choose>
    </xsl:for-each>
    free(element);
}

</xsl:if>
</xsl:template>

<xsl:template match="struct" mode="prototype">
<xsl:if test="@name='BlastOutput' or @name='Iteration'">
<xsl:value-of select="@name"/>Ptr parse<xsl:value-of select="@name"/>(xmlTextReaderPtr reader);
void dealloc<xsl:value-of select="@name"/>(<xsl:value-of select="@name"/>Ptr element);
</xsl:if> 
</xsl:template>

<xsl:template match="struct|field" mode="cname">
<xsl:value-of select="translate(@name,'-','_')"/>
</xsl:template>

</xsl:stylesheet>

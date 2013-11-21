from django import template
from django.template import Node
from django.template.base import TextNode

register = template.Library()

JS_CALL = '''
<script>
  $(function() {
    %(name)s(
      document.getElementById('JSELEMENT_%(id)s'),
      {%(subElements)s}
    );
  });
</script>
'''

ii = 0
def genID():
  global ii
  ii += 1
  return ii

class JSElementNode(Node):

  def __init__(self, name, id, subIDs, nodelist):
    self.name = name
    self.id = id
    self.subIDs = subIDs
    self.nodelist = nodelist

  def render(self, context):
    return self.nodelist.render(context) + self.jsCall(self.name)

  def jsCall(self, name):
    return JS_CALL % {
      'name': name,
      'id': self.id,
      'subElements': ','.join(
        '"%s": document.getElementById("JSELEMENT_%s")' % (k, v) for k, v in self.subIDs.iteritems()
      )
    }


@register.tag('jselem')
def jsElement(parser, token):
  bits = token.split_contents()
  if len(bits) < 2:
    raise Exception('{% jselem <ElementName> %} requires an element name!')
  name = bits[1]
  mainID = genID()
  subIDs = {}

  nodelist = parser.create_nodelist()
  while True:
    tmp_nodelist = parser.parse(('elemID', 'subID', 'endelem'))
    first_token = parser.tokens[0]
    token_contents = first_token.contents
    block_name = token_contents.split()[0]
    parser.delete_first_token()

    for node in tmp_nodelist:
      parser.extend_nodelist(nodelist, node, first_token)

    if block_name == 'elemID':
      node = TextNode('id="JSELEMENT_%d"' % mainID)
      parser.extend_nodelist(nodelist, node, first_token)
    elif block_name == 'subID':
      id = genID()
      subIDs[token_contents.split()[1]] = id
      node = TextNode('id="JSELEMENT_%d"' % id)
      parser.extend_nodelist(nodelist, node, first_token)
    else:
      break

  return JSElementNode(name, mainID, subIDs, nodelist)
